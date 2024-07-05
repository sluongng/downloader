use anyhow::{anyhow, Context, Result};
use clap::Parser;
use flate2::read::GzDecoder;
use serde::Deserialize;
use sha2::{Digest, Sha256};
use std::path::{Path, PathBuf};
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
}

#[derive(Deserialize)]
struct Config {
    urls: Vec<String>,
    digest: String,
    cache_prefix: String,
    file_type: Option<FileType>,
}

async fn load_config(config_path: &str) -> Result<Config> {
    let config_data = fs::read_to_string(config_path)
        .await
        .context("Unable to read config file")?;
    let config: Config =
        serde_json::from_str(&config_data).context("Unable to parse config file")?;
    Ok(config)
}

pub async fn main_with_args(args: Args) -> Result<()> {
    let config = load_config(&args.config).await?;
    let file_type = config
        .file_type
        .unwrap_or(file_type_from_urls(config.urls.clone()));

    download_from_urls(
        config.urls,
        &args.output,
        &config.digest,
        &config.cache_prefix,
        file_type,
    )
    .await?;
    println!("Download and validation successful!");

    Ok(())
}

#[derive(Deserialize)]
enum FileType {
    File,
    ArchiveZip,
    ArchiveTarGz,
}

fn file_type_from_urls(urls: Vec<String>) -> FileType {
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

async fn download_from_urls(
    urls: Vec<String>,
    output_path: &str,
    expected_digest: &str,
    cache_prefix: &str,
    file_type: FileType,
) -> Result<()> {
    // Check if the file is already cached
    let (hash_fn, hash, size) = digest_to_parts(expected_digest)?;
    let cache_path = Path::new(cache_prefix).join("cas").join(hash_fn).join(hash);
    if cache_path.exists() {
        // Copy the cached file or directory to the output
        copy_cache_to_output(&cache_path, output_path, file_type).await?;
        return Ok(());
    }

    for url in urls {
        match reqwest::get(&url).await {
            Ok(response) => {
                if response.status().is_success() {
                    let bytes = response.bytes().await?;

                    // Validate the digest
                    match validate_digest(&bytes, hash_fn, hash, size) {
                        Ok(_) => {
                            // Cache the file or directory
                            cache_bytes(&bytes, &cache_path, output_path, file_type).await?;
                            return Ok(());
                        }
                        Err(err) => {
                            println!("Digest validation failed for URL: {}. {}", url, err);
                        }
                    }
                } else {
                    println!(
                        "Failed to download from URL: {}. Status: {}",
                        url,
                        response.status()
                    );
                }
            }
            Err(e) => {
                println!("Error downloading from URL: {}. Error: {}", url, e);
            }
        }
    }
    Err(anyhow!("All URLs failed or digest validation failed"))
}

async fn copy_cache_to_output(
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
                        tokio::fs::set_permissions(&outpath, std::fs::Permissions::from_mode(mode))
                            .await?;
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

async fn cache_bytes(
    bytes: &[u8],
    cache_path: &PathBuf,
    output_path: &str,
    file_type: FileType,
) -> Result<()> {
    if let Some(parent) = cache_path.parent() {
        fs::create_dir_all(parent).await?;
    }
    let mut cache_file = File::create(&cache_path)
        .await
        .context("Failed to create cache file")?;
    cache_file
        .write_all(bytes)
        .await
        .context("Failed to write to cache file")?;
    cache_file
        .sync_all()
        .await
        .context("Failed to sync cache file")?;

    copy_cache_to_output(cache_path, output_path, file_type).await
}

fn digest_to_parts(digest: &str) -> Result<(&str, &str, usize)> {
    let (hash_fn, hash_and_size) = digest
        .split_once('-')
        .ok_or_else(|| anyhow!("Invalid digest format"))?;
    let (hash, size_str) = hash_and_size
        .split_once('/')
        .ok_or_else(|| anyhow!("Invalid digest format"))?;
    let size: usize = size_str
        .parse()
        .map_err(|_| anyhow!("Invalid digest size"))?;

    Ok((hash_fn, hash, size))
}

fn validate_digest(
    bytes: &[u8],
    hash_fn: &str,
    expected_hash: &str,
    expected_size: usize,
) -> Result<(), String> {
    let actual_size = bytes.len();
    let actual_hash = match hash_fn {
        "sha256" => format!("{:x}", Sha256::digest(bytes)),
        "blake3" => blake3::hash(bytes).to_hex().to_string(),
        _ => return Err(format!("Unsupported hash function: {}", hash_fn)),
    };

    if actual_hash == expected_hash && actual_size == expected_size {
        Ok(())
    } else {
        Err(format!(
            "Digest mismatch. Expected: {}/{}, Actual: {}/{}",
            expected_hash, expected_size, actual_hash, actual_size
        ))
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
        ($name:ident, $config:expr, $expected_value:expr) => {
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
                let digest = config["digest"].as_str().unwrap();
                let (hash_fn, hash, _) = digest_to_parts(digest)?;
                let cached_file_path = cache_path.join("cas").join(hash_fn).join(hash);

                // Run the download
                let args = Args {
                    config: config_path.to_str().unwrap().to_string(),
                    output: output_path.to_str().unwrap().to_string(),
                };
                main_with_args(args).await?;

                assert!(fs::metadata(&output_path).is_ok());
                let output_content = fs::read_to_string(&output_path)?;
                assert_eq!(output_content, $expected_value);

                assert!(fs::metadata(&cached_file_path).is_ok());
                let cached_content = fs::read_to_string(&cached_file_path)?;
                assert_eq!(cached_content, $expected_value);

                // Send the shutdown signal to the server
                let _ = shutdown_tx.send(());
                fs::remove_file(&output_path)?;

                // Run the download again to test for cache
                let args = Args {
                    config: config_path.to_str().unwrap().to_string(),
                    output: output_path.to_str().unwrap().to_string(),
                };
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
        ($name:ident, $config:expr, $expected_structure:expr) => {
            #[tokio::test]
            async fn $name() -> Result<()> {
                let config = $config;
                let expected_structure = $expected_structure;

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
                let digest = config["digest"].as_str().unwrap();
                let (hash_fn, hash, _) = digest_to_parts(digest)?;
                let cached_file_path = cache_path.join("cas").join(hash_fn).join(hash);

                // Run the download
                let args = Args {
                    config: config_path.to_str().unwrap().to_string(),
                    output: output_path.to_str().unwrap().to_string(),
                };
                main_with_args(args).await?;

                // Verify the extracted structure
                verify_directory_structure(&output_path, &expected_structure)?;

                // Check the cache directory contents
                assert!(cached_file_path.is_file(), "Cached file does not exist");

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

    test_download!(
        base_case,
        json!({
            "urls": ["http://localhost:{port}/file_to_download.txt"],
            "digest": "sha256-315f5bdb76d078c43b8ac0064e4a0164612b1fce77c869345bfc94c75894edd3/13",
            "cache_prefix": "{cache_prefix}"
        }),
        "Hello, world!"
    );
    test_download!(
        with_different_content,
        json!({
            "urls": ["http://localhost:{port}/file_to_download.txt"],
            "digest": "sha256-9ca7a7d04afda0e4e47db2dbdd9878903674918451e7963fb42777a6a84bdd4a/18",
            "cache_prefix": "{cache_prefix}"
        }),
        "Hello, downloader!"
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
        "Hello, world!"
    );
    test_download!(
        with_blake3,
        json!({
            "urls": ["http://localhost:{port}/file_to_download.txt"],
            "digest": "blake3-ede5c0b10f2ec4979c69b52f61e42ff5b413519ce09be0f14d098dcfe5f6f98d/13",
            "cache_prefix": "{cache_prefix}"
        }),
        "Hello, world!"
    );

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
        ]
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
        ]
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
        ]
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
        ($name:ident, $config:expr, $expected_structure:expr) => {
            #[tokio::test]
            async fn $name() -> Result<()> {
                let config = $config;
                let expected_structure = $expected_structure;

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
                let digest = config["digest"].as_str().unwrap();
                let (hash_fn, hash, _) = digest_to_parts(digest)?;
                let cached_file_path = cache_path.join("cas").join(hash_fn).join(hash);

                // Run the download
                let args = Args {
                    config: config_path.to_str().unwrap().to_string(),
                    output: output_path.to_str().unwrap().to_string(),
                };
                main_with_args(args).await?;

                // Verify the extracted structure
                verify_directory_structure(&output_path, &expected_structure)?;

                // Check the cache directory contents
                assert!(cached_file_path.is_file(), "Cached file does not exist");

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
        ]
    );
}
