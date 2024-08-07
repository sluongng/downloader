use anyhow::{anyhow, Context, Result};
use clap::Parser;
use core::fmt;
use flate2::read::GzDecoder;
use serde::Deserialize;
use sha2::{Digest, Sha256};
use std::path::{Path, PathBuf};
use std::str::FromStr;
use tar::Archive;
use tokio::fs::{self, File};
use tokio::io::AsyncWriteExt;
use zip::ZipArchive;

#[derive(Parser, Debug, Clone)]
#[clap(author, version, about, long_about = None)]
pub struct Args {
    /// Path to the JSON config file
    #[arg(short, long)]
    pub config: String,

    /// Output file name or directory
    #[arg(short, long)]
    pub output: String,

    /// Cache directory (optional)
    #[arg(short, long)]
    pub cache: Option<PathBuf>,
}

#[derive(Deserialize)]
struct Config {
    urls: Vec<String>,
    digest: String,
}

struct Downloader {
    cache: Cache,
    verifier: Verifier,
}

impl Downloader {
    fn new(cache_base: Option<PathBuf>, verifier: Verifier) -> Self {
        let cache_base_path =
            cache_base.unwrap_or_else(|| dirs::cache_dir().expect("could not detect cache dir"));

        Self {
            cache: Cache::new(cache_base_path),
            verifier,
        }
    }

    fn cache_path(&self) -> PathBuf {
        self.cache.get_path(
            self.verifier.digest_fn.to_string().as_str(),
            &self.verifier.hash,
        )
    }

    async fn download_from_urls(&self, urls: Vec<String>, output_path: &str) -> Result<()> {
        let file_type = FileType::from(urls.clone());
        let cache_path = self.cache_path();
        if cache_path.exists() {
            self.cache
                .copy_to_output(&cache_path, output_path, file_type)
                .await?;
            return Ok(());
        }

        for url in urls {
            match reqwest::get(&url).await.context("Failed to get URL") {
                Ok(response) => {
                    if response.status().is_success() {
                        let bytes = response.bytes().await?;
                        if self.verifier.validate_digest(&bytes).is_ok() {
                            self.cache.write(&cache_path, &bytes).await?;
                            return self
                                .cache
                                .copy_to_output(&cache_path, output_path, file_type)
                                .await
                                .context("Failed to copy to output after successful download");
                        }
                    }
                }
                Err(e) => println!("Error downloading from URL: {}. Error: {}", url, e),
            }
        }
        Err(anyhow!("All URLs failed or digest validation failed"))
    }
}

struct Cache {
    base_path: PathBuf,
}

impl Cache {
    fn new(base_path: PathBuf) -> Self {
        Self { base_path }
    }

    fn get_path(&self, hash_fn: &str, hash: &str) -> PathBuf {
        self.base_path
            .join("downloader")
            .join("cas")
            .join(hash_fn)
            .join(hash)
    }

    async fn copy_to_output(
        &self,
        cache_path: &PathBuf,
        output_path: &str,
        file_type: FileType,
    ) -> Result<()> {
        match file_type {
            FileType::ArchiveZip => {
                fs::create_dir_all(output_path).await?;
                let cache_file = File::open(&cache_path)
                    .await
                    .context("Failed to open cache file")?;
                let mut archive = ZipArchive::new(cache_file.into_std().await)?;

                tokio::fs::create_dir_all(output_path).await?;
                for i in 0..archive.len() {
                    let mut file = archive.by_index(i)?;
                    let outpath = Path::new(output_path).join(file.mangled_name());

                    if file.name().ends_with('/') {
                        tokio::fs::create_dir_all(&outpath).await?;
                    } else {
                        if let Some(p) = outpath.parent() {
                            if !p.exists() {
                                tokio::fs::create_dir_all(p).await?;
                            }
                        }
                        let outfile = File::create(&outpath).await?;
                        std::io::copy(&mut file, &mut outfile.into_std().await)?;
                    }

                    // Get and Set permissions
                    #[cfg(unix)]
                    {
                        use std::os::unix::fs::PermissionsExt;
                        if let Some(mode) = file.unix_mode() {
                            tokio::fs::set_permissions(
                                &outpath,
                                std::fs::Permissions::from_mode(mode),
                            )
                            .await
                            .context("failed to set file permission")?;
                        }
                    }
                }
            }
            FileType::ArchiveTarGz => {
                let cache_file = File::open(&cache_path)
                    .await
                    .context("Failed to open cache file")?;
                let cache_file_std = cache_file.into_std().await;
                let tar = GzDecoder::new(cache_file_std);
                let mut archive = Archive::new(tar);

                for entry in archive.entries()? {
                    let mut entry = entry?;
                    let path = entry.path()?;
                    let outpath = Path::new(output_path).join(path);

                    if let Some(p) = outpath.parent() {
                        if !p.exists() {
                            tokio::fs::create_dir_all(p).await?;
                        }
                    }

                    entry.unpack(&outpath)?;
                }
            }
            FileType::File => {
                let mut cache_file = File::open(&cache_path)
                    .await
                    .context("Failed to open cache file")?;
                let mut output_file = File::create(output_path)
                    .await
                    .context("Failed to create output file")?;
                tokio::io::copy(&mut cache_file, &mut output_file)
                    .await
                    .context("Failed to copy from cache to output")?;
                output_file
                    .sync_all()
                    .await
                    .context("Failed to sync output file")?;
            }
        }
        Ok(())
    }

    async fn write(&self, cache_path: &PathBuf, bytes: &[u8]) -> Result<()> {
        if let Some(parent) = cache_path.parent() {
            fs::create_dir_all(parent).await?;
        }
        println!("writing cache file: {:?}", cache_path);
        let mut cache_file = File::create(&cache_path)
            .await
            .context("Failed to create cache file")?;
        cache_file
            .write_all(bytes)
            .await
            .context("Failed to write to cache file")?;
        cache_file
            .flush()
            .await
            .context("Failed to sync cache file")
    }
}

#[derive(Debug, PartialEq)]
enum DigestFn {
    SHA256,
    BLAKE3,
}

#[derive(Debug, PartialEq, Eq)]
struct ParseDigestFnError;

impl FromStr for DigestFn {
    type Err = ParseDigestFnError;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "sha256" => Ok(DigestFn::SHA256),
            "blake3" => Ok(DigestFn::BLAKE3),
            _ => Err(ParseDigestFnError),
        }
    }
}

impl fmt::Display for DigestFn {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", format!("{:?}", self).to_lowercase())
    }
}

impl DigestFn {
    fn compute_digest(&self, bytes: &[u8]) -> (String, usize) {
        let size = bytes.len();
        let hash = match self {
            DigestFn::SHA256 => format!("{:x}", Sha256::digest(bytes)),
            DigestFn::BLAKE3 => blake3::hash(bytes).to_hex().to_string(),
        };
        (hash, size)
    }
}

struct Verifier {
    digest_fn: DigestFn,
    hash: String,
    size: usize,
}

impl FromStr for Verifier {
    type Err = ();

    fn from_str(digest: &str) -> Result<Self, Self::Err> {
        let (hash_fn, hash_and_size) = digest
            .split_once('-')
            .ok_or_else(|| anyhow!("Invalid digest format"))
            .unwrap();
        let digest_fn = DigestFn::from_str(hash_fn)
            .unwrap_or_else(|_| panic!("Invalid digest function {}", hash_fn));
        let (hash, size_str) = hash_and_size
            .split_once('/')
            .ok_or_else(|| anyhow!("Invalid digest format"))
            .unwrap();
        let size: usize = size_str
            .parse()
            .map_err(|_| anyhow!("Invalid digest size"))
            .unwrap();
        let hash = hash.to_string();

        Ok(Self {
            digest_fn,
            hash,
            size,
        })
    }
}

impl Verifier {
    fn validate_digest(&self, bytes: &[u8]) -> Result<(), String> {
        let (actual_hash, actual_size) = self.digest_fn.compute_digest(bytes);

        if actual_hash == self.hash && actual_size == self.size {
            Ok(())
        } else {
            Err(format!(
                "Digest mismatch. Expected: {}/{}, Actual: {}/{}",
                self.hash, self.size, actual_hash, actual_size
            ))
        }
    }
}

async fn load_config(config_path: &str) -> Result<Config> {
    let config_data = fs::read_to_string(config_path)
        .await
        .context("Unable to read config file")?;
    serde_json::from_str(&config_data).context("Unable to parse config file")
}

pub async fn main_with_args(args: Args) -> Result<()> {
    let config = load_config(&args.config).await?;
    let verifier = Verifier::from_str(&config.digest).unwrap();
    Downloader::new(args.cache, verifier)
        .download_from_urls(config.urls, &args.output)
        .await
}

#[derive(Deserialize)]
enum FileType {
    File,
    ArchiveZip,
    ArchiveTarGz,
}

impl From<Vec<String>> for FileType {
    fn from(urls: Vec<String>) -> Self {
        for url in urls {
            if url.ends_with(".zip") {
                return FileType::ArchiveZip;
            }
            if url.ends_with("tar.gz") {
                return FileType::ArchiveTarGz;
            }
        }
        FileType::File
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use anyhow::Result;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use serde_json::json;
    use std::fs::{self, File};
    use std::io::Cursor;
    use std::io::Write;
    use tar::Builder;
    use tempfile::tempdir;
    use tokio::sync::oneshot;
    use warp::Filter;
    use zip::write::SimpleFileOptions;
    use zip::{DateTime, ZipWriter};

    async fn start_server(expected_value: &'static str) -> (u16, oneshot::Sender<()>) {
        let route = warp::path("file_to_download.txt")
            .map(move || warp::reply::with_header(expected_value, "content-type", "text/plain"));

        let (shutdown_tx, shutdown_rx) = oneshot::channel();
        let (addr, server) =
            warp::serve(route).bind_with_graceful_shutdown(([127, 0, 0, 1], 0), async {
                shutdown_rx.await.ok();
            });

        tokio::spawn(server);
        (addr.port(), shutdown_tx)
    }

    macro_rules! test_download {
        ($name:ident, $config:expr, $expected_value:expr, $expected_cache_path:expr) => {
            #[tokio::test]
            async fn $name() -> Result<()> {
                let config = $config;

                // Start the local server
                let (port, shutdown_tx) = start_server($expected_value).await;

                // Set up a temporary directory for test environment
                let tmp_dir = tempdir()?;
                let output_path = tmp_dir.path().join("output.txt");
                let cache_path = tmp_dir.path().join("cache");

                // Create config file
                let config_path = tmp_dir.path().join("config.json");
                let mut config_file = File::create(&config_path)?;
                let config_str = config
                    .to_string()
                    .replace("{port}", &port.to_string())
                    .replace("{cache_prefix}", cache_path.to_str().unwrap());
                write!(config_file, "{}", config_str)?;

                // Extract the hash from the digest
                let cached_file_path = $expected_cache_path
                    .to_string()
                    .replace("{cache_prefix}", cache_path.to_str().unwrap());

                // Run the download
                let args = Args {
                    config: config_path.to_str().unwrap().to_string(),
                    output: output_path.to_str().unwrap().to_string(),
                    cache: Some(PathBuf::from(cache_path)),
                };
                main_with_args(args.clone()).await?;

                assert!(fs::metadata(&output_path).is_ok());
                let output_content = fs::read_to_string(&output_path)?;
                assert_eq!(output_content, $expected_value);

                assert!(
                    fs::metadata(&cached_file_path).is_ok(),
                    "Cached file {:?} does not exist",
                    cached_file_path
                );
                let cached_content = fs::read_to_string(&cached_file_path)?;
                assert_eq!(cached_content, $expected_value);

                // Send the shutdown signal to the server
                let _ = shutdown_tx.send(());
                fs::remove_file(&output_path)?;

                // Run the download again to test for cache
                main_with_args(args).await?;

                assert!(fs::metadata(&output_path).is_ok());
                let output_content = fs::read_to_string(&output_path)?;
                assert_eq!(output_content, $expected_value);

                assert!(fs::metadata(&cached_file_path).is_ok());
                let cached_content = fs::read_to_string(&cached_file_path)?;
                assert_eq!(cached_content, $expected_value);

                Ok(())
            }
        };
    }
    test_download!(
        base_case,
        json!({
            "urls": ["http://localhost:{port}/file_to_download.txt"],
            "digest": "sha256-315f5bdb76d078c43b8ac0064e4a0164612b1fce77c869345bfc94c75894edd3/13",
            "cache_prefix": "{cache_prefix}"
        }),
        "Hello, world!",
        "{cache_prefix}/downloader/cas/sha256/315f5bdb76d078c43b8ac0064e4a0164612b1fce77c869345bfc94c75894edd3"
    );
    test_download!(
        with_different_content,
        json!({
            "urls": ["http://localhost:{port}/file_to_download.txt"],
            "digest": "sha256-9ca7a7d04afda0e4e47db2dbdd9878903674918451e7963fb42777a6a84bdd4a/18",
            "cache_prefix": "{cache_prefix}"
        }),
        "Hello, downloader!",
        "{cache_prefix}/downloader/cas/sha256/9ca7a7d04afda0e4e47db2dbdd9878903674918451e7963fb42777a6a84bdd4a"
    );
    test_download!(
        with_invalid_url,
        json!({
            "urls": [
                "http://invalid.url/somefile.txt",
                "http://localhost:{port}/file_to_download.txt",
            ],
            "digest": "sha256-315f5bdb76d078c43b8ac0064e4a0164612b1fce77c869345bfc94c75894edd3/13",
            "cache_prefix": "{cache_prefix}"
        }),
        "Hello, world!",
        "{cache_prefix}/downloader/cas/sha256/315f5bdb76d078c43b8ac0064e4a0164612b1fce77c869345bfc94c75894edd3"
    );
    test_download!(
        with_blake3,
        json!({
            "urls": ["http://localhost:{port}/file_to_download.txt"],
            "digest": "blake3-ede5c0b10f2ec4979c69b52f61e42ff5b413519ce09be0f14d098dcfe5f6f98d/13",
            "cache_prefix": "{cache_prefix}"
        }),
        "Hello, world!",
        "{cache_prefix}/downloader/cas/blake3/ede5c0b10f2ec4979c69b52f61e42ff5b413519ce09be0f14d098dcfe5f6f98d"
    );

    // Helper function to start a server with custom zip content
    async fn start_custom_zip_server(
        structure: &[ArchiveEntry],
    ) -> Result<(u16, oneshot::Sender<()>)> {
        let zip_data = generate_custom_zip_data(structure)?;
        let route = warp::path("archive_to_download.zip").map(move || {
            warp::reply::with_header(zip_data.clone(), "content-type", "application/zip")
        });

        let (shutdown_tx, shutdown_rx) = oneshot::channel();
        let (addr, server) =
            warp::serve(route).bind_with_graceful_shutdown(([127, 0, 0, 1], 0), async {
                shutdown_rx.await.ok();
            });

        tokio::spawn(server);
        Ok((addr.port(), shutdown_tx))
    }

    fn generate_custom_zip_data(structure: &[ArchiveEntry]) -> Result<Vec<u8>> {
        let mut buffer = Vec::new();
        let mut zip = ZipWriter::new(Cursor::new(&mut buffer));
        let options = SimpleFileOptions::default()
            .compression_method(zip::CompressionMethod::Stored)
            .last_modified_time(DateTime::from_date_and_time(1980, 1, 1, 0, 0, 0)?)
            .unix_permissions(0o755);

        fn add_entries(
            zip: &mut ZipWriter<Cursor<&mut Vec<u8>>>,
            structure: &[ArchiveEntry],
            path: &str,
            options: &SimpleFileOptions,
        ) -> Result<()> {
            for entry in structure {
                match entry {
                    ArchiveEntry::File(name) => {
                        let full_path = format!("{}{}", path, name);
                        zip.start_file(full_path.clone(), *options)?;
                        zip.write_all(format!("Content of {}", name).as_bytes())?;
                    }
                    ArchiveEntry::Directory(name, children) => {
                        let full_path = format!("{}{}/", path, name);
                        zip.add_directory(&full_path, *options)?;
                        add_entries(zip, children, &full_path, options)?;
                    }
                }
            }
            Ok(())
        }

        add_entries(&mut zip, structure, "", &options)?;
        zip.finish()?;
        Ok(buffer)
    }

    #[derive(Debug, Clone)]
    enum ArchiveEntry {
        File(String),
        Directory(String, Vec<ArchiveEntry>),
    }

    fn count_entries(structure: &[ArchiveEntry]) -> usize {
        structure.iter().fold(0, |acc, entry| match entry {
            ArchiveEntry::File(_) => acc + 1,
            ArchiveEntry::Directory(_, children) => acc + 1 + count_entries(children),
        })
    }

    fn verify_directory_structure(base_path: &Path, expected: &[ArchiveEntry]) -> Result<()> {
        for entry in expected {
            match entry {
                ArchiveEntry::File(name) => {
                    let file_path = base_path.join(name);
                    assert!(
                        file_path.is_file(),
                        "Expected file does not exist: {:?}",
                        file_path
                    );
                    let content = fs::read_to_string(&file_path)?;
                    assert_eq!(
                        content,
                        format!("Content of {}", name),
                        "File content mismatch for {:?}",
                        file_path
                    );
                }
                ArchiveEntry::Directory(name, children) => {
                    let dir_path = base_path.join(name);
                    assert!(
                        dir_path.is_dir(),
                        "Expected directory does not exist: {:?}",
                        dir_path
                    );
                    verify_directory_structure(&dir_path, children)?;
                }
            }
        }
        Ok(())
    }

    macro_rules! test_download_archive {
        ($name:ident, $config:expr, $expected_content:expr, $expected_cache_path:expr) => {
            #[tokio::test]
            async fn $name() -> Result<()> {
                let config = $config;
                let expected_structure = $expected_content;

                // Start the local server with custom zip content
                let (port, shutdown_tx) = start_custom_zip_server(&expected_structure).await?;

                // Set up a temporary directory for test environment
                let tmp_dir = tempdir()?;
                let output_path = tmp_dir.path().join("output_dir");
                let cache_path = tmp_dir.path().join("cache");

                // Create config file
                let config_path = tmp_dir.path().join("config.json");
                let mut config_file = File::create(&config_path)?;
                let config_str = config
                    .to_string()
                    .replace("{port}", &port.to_string())
                    .replace("{cache_prefix}", cache_path.to_str().unwrap());
                write!(config_file, "{}", config_str)?;

                // Extract the hash from the digest
                let cached_file_path = $expected_cache_path
                    .to_string()
                    .replace("{cache_prefix}", cache_path.to_str().unwrap());

                // Run the download
                let args = Args {
                    config: config_path.to_str().unwrap().to_string(),
                    output: output_path.to_str().unwrap().to_string(),
                    cache: Some(PathBuf::from(cache_path)),
                };
                main_with_args(args).await?;

                // Verify the extracted structure
                verify_directory_structure(&output_path, &expected_structure)?;

                // Check the cache archive file exist
                assert!(
                    fs::metadata(&cached_file_path).is_ok(),
                    "Cached file {:?} does not exist",
                    cached_file_path
                );

                // Verify the cached file is a valid zip
                let cached_file = File::open(&cached_file_path)?;
                let archive = ZipArchive::new(cached_file)?;
                let expected_entry_count = count_entries(&expected_structure);
                assert_eq!(
                    archive.len(),
                    expected_entry_count,
                    "Cached zip file should contain {} entries",
                    expected_entry_count
                );

                // Send the shutdown signal to the server
                let _ = shutdown_tx.send(());

                Ok(())
            }
        };
    }

    test_download_archive!(
        zip_with_files,
        json!({
            "urls": ["http://localhost:{port}/archive_to_download.zip"],
            "digest": "sha256-f3d4199459e64d5b8c2e5ebc426970a1c1cc4b3b0d704789fb0975f4cc180644/250",
            "cache_prefix": "{cache_prefix}"
        }),
        vec![
            ArchiveEntry::File("file1.txt".to_string()),
            ArchiveEntry::File("file2.txt".to_string())
        ],
        "{cache_prefix}/downloader/cas/sha256/f3d4199459e64d5b8c2e5ebc426970a1c1cc4b3b0d704789fb0975f4cc180644"
    );
    test_download_archive!(
        zip_with_dir,
        json!({
            "urls": ["http://localhost:{port}/archive_to_download.zip"],
            "digest": "sha256-dedd5459a6b30a8f74836d562a13398a96d4d1e56443aa8de74069c6b2ec5ac1/342",
            "cache_prefix": "{cache_prefix}"
        }),
        vec![
            ArchiveEntry::File("file1.txt".to_string()),
            ArchiveEntry::Directory(
                "foo".to_string(),
                vec![ArchiveEntry::File("file2.txt".to_string())]
            ),
        ],
        "{cache_prefix}/downloader/cas/sha256/dedd5459a6b30a8f74836d562a13398a96d4d1e56443aa8de74069c6b2ec5ac1"
    );
    test_download_archive!(
        zip_with_blake3,
        json!({
            "urls": ["http://localhost:{port}/archive_to_download.zip"],
            "digest": "blake3-c03560f38064a39ae621b84e898b689f9d41bbf6cf0af49895ff2bf5f626d3e9/342",
            "cache_prefix": "{cache_prefix}"
        }),
        vec![
            ArchiveEntry::File("file1.txt".to_string()),
            ArchiveEntry::Directory(
                "foo".to_string(),
                vec![ArchiveEntry::File("file2.txt".to_string())]
            ),
        ],
        "{cache_prefix}/downloader/cas/blake3/c03560f38064a39ae621b84e898b689f9d41bbf6cf0af49895ff2bf5f626d3e9"
    );

    // Helper function to start a server with custom tar.gz content
    async fn start_custom_targz_server(
        structure: &[ArchiveEntry],
    ) -> Result<(u16, oneshot::Sender<()>)> {
        let tar_gz_data = generate_custom_targz_data(structure)?;
        let route = warp::path("archive_to_download.tar.gz").map(move || {
            warp::reply::with_header(tar_gz_data.clone(), "content-type", "application/gzip")
        });

        let (shutdown_tx, shutdown_rx) = oneshot::channel();
        let (addr, server) =
            warp::serve(route).bind_with_graceful_shutdown(([127, 0, 0, 1], 0), async {
                shutdown_rx.await.ok();
            });

        tokio::spawn(server);
        Ok((addr.port(), shutdown_tx))
    }

    fn generate_custom_targz_data(structure: &[ArchiveEntry]) -> Result<Vec<u8>> {
        let mut buffer = Vec::new();
        let encoder = GzEncoder::new(&mut buffer, Compression::default());
        let mut tar = Builder::new(encoder);

        fn add_entries(
            tar: &mut Builder<GzEncoder<&mut Vec<u8>>>,
            structure: &[ArchiveEntry],
            path: &str,
        ) -> Result<()> {
            for entry in structure {
                match entry {
                    ArchiveEntry::File(name) => {
                        let full_path = format!("{}{}", path, name);
                        let data = format!("Content of {}", name);
                        let mut header = tar::Header::new_old();
                        header.set_size(data.len() as u64);
                        tar.append_data(&mut header, full_path, data.as_bytes())?;
                    }
                    ArchiveEntry::Directory(name, children) => {
                        let full_path = format!("{}{}/", path, name);
                        tar.append_dir_all(full_path.clone(), full_path.clone())?;
                        add_entries(tar, children, &format!("{}{}", path, name))?;
                    }
                }
            }
            Ok(())
        }

        add_entries(&mut tar, structure, "")?;
        tar.into_inner()?.finish()?;
        Ok(buffer)
    }

    macro_rules! test_download_targz {
        ($name:ident, $config:expr, $expected_content:expr, $expected_cache_path:expr) => {
            #[tokio::test]
            async fn $name() -> Result<()> {
                let config = $config;
                let expected_structure = $expected_content;

                // Start the local server with custom tar.gz content
                let (port, shutdown_tx) = start_custom_targz_server(&expected_structure).await?;

                // Set up a temporary directory for test environment
                let tmp_dir = tempdir()?;
                let output_path = tmp_dir.path().join("output_dir");
                let cache_path = tmp_dir.path().join("cache");

                // Create config file
                let config_path = tmp_dir.path().join("config.json");
                let mut config_file = File::create(&config_path)?;
                let config_str = config
                    .to_string()
                    .replace("{port}", &port.to_string())
                    .replace("{cache_prefix}", cache_path.to_str().unwrap());
                write!(config_file, "{}", config_str)?;

                // Extract the hash from the digest
                let cached_file_path = $expected_cache_path
                    .to_string()
                    .replace("{cache_prefix}", cache_path.to_str().unwrap());

                // Run the download
                let args = Args {
                    config: config_path.to_str().unwrap().to_string(),
                    output: output_path.to_str().unwrap().to_string(),
                    cache: Some(PathBuf::from(cache_path)),
                };
                main_with_args(args).await?;

                // Verify the extracted structure
                verify_directory_structure(&output_path, &expected_structure)?;

                // Check the cached archive file contents
                assert!(
                    fs::metadata(&cached_file_path).is_ok(),
                    "Cached file {:?} does not exist",
                    cached_file_path
                );

                // Verify the cached file is a valid tar.gz
                let cached_file = File::open(&cached_file_path)?;
                let tar_gz = GzDecoder::new(cached_file);
                let mut archive = Archive::new(tar_gz);
                let expected_entry_count = count_entries(&expected_structure);

                assert_eq!(
                    archive.entries()?.count(),
                    expected_entry_count,
                    "Cached tar.gz file should contain {} entries",
                    expected_entry_count
                );

                // Send the shutdown signal to the server
                let _ = shutdown_tx.send(());

                Ok(())
            }
        };
    }
    test_download_targz!(
        targz_with_files,
        json!({
            "urls": ["http://localhost:{port}/archive_to_download.tar.gz"],
            "digest": "sha256-51d2ee272e361eecafef092b4330df7c40ebac3375212a1b5d663217d9f3d4a4/101",
            "cache_prefix": "{cache_prefix}"
        }),
        vec![
            ArchiveEntry::File("file1.txt".to_string()),
            ArchiveEntry::File("file2.txt".to_string())
        ],
        "{cache_prefix}/downloader/cas/sha256/51d2ee272e361eecafef092b4330df7c40ebac3375212a1b5d663217d9f3d4a4"
    );
}
