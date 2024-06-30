use anyhow::Result;
use clap::Parser;
use downloader::Args;

#[tokio::main]
async fn main() -> Result<()> {
    let args = Args::parse();
    downloader::main_with_args(args).await
}
