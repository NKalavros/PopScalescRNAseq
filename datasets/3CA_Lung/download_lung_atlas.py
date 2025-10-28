#!/usr/bin/env python3
"""Download Curated Cancer Cell Atlas - Lung Cancer datasets.

This script downloads single-cell RNA-seq data and metadata from the
Curated Cancer Cell Atlas (3CA) Lung Cancer collection hosted at Weizmann Institute.

Source: https://www.weizmann.ac.il/sites/3CA/lung

Features
========
* Download individual studies or all studies
* Download data, metadata, or both
* Resume capability for interrupted downloads
* Progress tracking
* Automatic extraction of tar.gz files
* Verify downloads with file size checks

Usage
=====
Download all data and metadata:
    python download_lung_atlas.py --output-dir ./LungAtlas --download-all

Download specific studies:
    python download_lung_atlas.py --output-dir ./LungAtlas --studies Bischoff2021 Maynard2020

Download only metadata for all studies:
    python download_lung_atlas.py --output-dir ./LungAtlas --download-all --metadata-only

List available studies:
    python download_lung_atlas.py --list-studies
"""

from __future__ import annotations

import argparse
import hashlib
import logging
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional
from urllib.request import urlopen, Request
from urllib.error import URLError, HTTPError

logger = logging.getLogger(__name__)

# Dataset definitions
DATASETS = {
    "Bischoff2021": {
        "full_name": "Bischoff et al. 2021",
        "description": "Lung adenocarcinoma (10x, 20 samples, 120,961 cells)",
        "data_url": "https://www.dropbox.com/scl/fi/08bjhmr5b7zcr9w3qm37n/Data_Bischoff2021_Lung.tar.gz?rlkey=xlv2vwb5xenzczoo9jiktr1ry&dl=1",
        "metadata_url": "https://www.dropbox.com/scl/fi/5phcsqxntjdjjfxypf4l3/Meta-data_Bischoff2021_Lung.tar.gz?rlkey=gmn535v9fqtqaqibvjria1lks&dl=1",
    },
    "Chan2021": {
        "full_name": "Chan et al. 2021",
        "description": "SCLC/Lung Adenocarcinoma (10x, 37 samples, 86,662 cells)",
        "data_url": "https://www.dropbox.com/scl/fi/mn5xe6x5hufjy2j570r3y/Data_Chan2021_Lung.tar.gz?rlkey=gptwkhhao0yz3xw0b0olf1l9u&dl=1",
        "metadata_url": "https://www.dropbox.com/scl/fi/n6i0d3km34y5txol8hicy/Meta-data_Chan2021_Lung.tar.gz?rlkey=c6wsubisf8jsnq4ydsppkueec&dl=1",
    },
    "Guo2018": {
        "full_name": "Guo et al. 2018",
        "description": "NSCLC (SmartSeq2, 40 samples, 12,346 cells)",
        "data_url": "https://www.dropbox.com/scl/fi/xm4x3thusfqi7vhqbkds6/Data_Guo2018_Lung.tar.gz?rlkey=cnxhnn4c7brjbq6cbdbo2cg5u&dl=1",
        "metadata_url": "https://www.dropbox.com/scl/fi/ni6ja6vcahj1xehyesy6w/Meta-data_Guo2018_Lung.tar.gz?rlkey=dsyaeztbjoa4zkz5p6vjp1ya1&dl=1",
    },
    "Ireland2020": {
        "full_name": "Ireland et al. 2020",
        "description": "SCLC (10x, 12 samples, 31,258 cells)",
        "data_url": "https://www.dropbox.com/scl/fi/8b4p8d551dv9wp0mt4blg/Data_Ireland2020_Lung.tar.gz?rlkey=xidxajfu4ynnbl6f3lw9owqlm&dl=1",
        "metadata_url": "https://www.dropbox.com/scl/fi/0k2lt05cdvc6639h2j3uq/Meta-data_Ireland2020_Lung.tar.gz?rlkey=t9jwymniruu2ia4fh0s33n4di&dl=1",
    },
    "Kim2020": {
        "full_name": "Kim et al. 2020",
        "description": "Lung adenocarcinoma (10x, 14 samples, 32,493 cells)",
        "data_url": "https://www.dropbox.com/scl/fi/n3hg4fj9bth2tkroej6t2/Data_Kim2020_Lung.tar.gz?rlkey=fxabtlt57pexdj7hdhvddzyjx&dl=1",
        "metadata_url": "https://www.dropbox.com/scl/fi/0vvy8pugcvfz50cut60gy/Meta-data_Kim2020_Lung.tar.gz?rlkey=hrg78tcmlt7l7py1em19y51g7&dl=1",
    },
    "Laughney2020": {
        "full_name": "Laughney et al. 2020",
        "description": "Lung cancer (10x, 17 samples, 40,505 cells)",
        "data_url": "https://www.dropbox.com/scl/fi/idlzf638lbwtejs3oruoe/Data_Laughney2020_Lung.tar.gz?rlkey=wu1a2fbj2n4iyya5mf2q1on1b&dl=1",
        "metadata_url": "https://www.dropbox.com/scl/fi/bxagrx8emjcca2h69qwk0/Meta-data_Laughney2020_Lung.tar.gz?rlkey=b715kxnv7f7edc11itp0ap417&dl=1",
    },
    "Maynard2020": {
        "full_name": "Maynard et al. 2020",
        "description": "Lung cancer (SmartSeq2, 50 samples, 19,777 cells)",
        "data_url": "https://www.dropbox.com/scl/fi/py7z9tgjdp1m3aw2zop7x/Data_Maynard2020_Lung.tar.gz?rlkey=rirsnfhqus3indqpp5g02ce1u&dl=1",
        "metadata_url": "https://www.dropbox.com/scl/fi/r62e9vf0uo60i5fpfsrd3/Meta-data_Maynard2020_Lung.tar.gz?rlkey=on1mc7puyky0tsnlzkg12l0zb&dl=1",
    },
    "Qian2020": {
        "full_name": "Qian et al. 2020",
        "description": "Lung cancer (10x, 8 samples, 27,262 cells)",
        "data_url": "https://www.dropbox.com/scl/fi/361w3739nwvc89l2gx4u8/Data_Qian2020_Lung.tar.gz?rlkey=4w00kvwlnf7hdhvddzyjx&dl=1",
        "metadata_url": "https://www.dropbox.com/scl/fi/4jtax8aklbhszx80krl3f/Meta-data_Qian2020_Lung.tar.gz?rlkey=6hbdtx9ro7n8hdi0z5jti0ikd&dl=1",
    },
    "Song2019": {
        "full_name": "Song et al. 2019",
        "description": "NSCLC (10x, 6 samples, 8,772 cells)",
        "data_url": "https://www.dropbox.com/scl/fi/k7yd3cine5dmwsjgv839h/Data_Song2019_Lung.tar.gz?rlkey=y9ipfcxsu5btza7j8n2wvpnul&dl=1",
        "metadata_url": "https://www.dropbox.com/scl/fi/eyecfgjefrksfv3ve0qar/Meta-data_Song2019_Lung.tar.gz?rlkey=qiz9mluost6nqh11ddcze84w2&dl=1",
    },
    "Xing2021": {
        "full_name": "Xing et al. 2021",
        "description": "Lung cancer (10x, 31 samples, 118,293 cells)",
        "data_url": "https://www.dropbox.com/scl/fi/3ll4j9tdur4485of8s3yf/Data_Xing2021_Lung.tar.gz?rlkey=a4l934jl6hwgqo2tx8f4fuxrk&dl=1",
        "metadata_url": "https://www.dropbox.com/scl/fi/7rawj2jnq8m4wthh9alvg/Meta-data_Xing2021_Lung.tar.gz?rlkey=fdj0tbhi1zi8gnkf6ng724m8t&dl=1",
    },
    "Zilionis2019": {
        "full_name": "Zilionis et al. 2019",
        "description": "NSCLC (inDrop, 7 samples, 31,179 cells)",
        "data_url": "https://www.dropbox.com/scl/fi/ot0tjg6254ne6hqzp11zv/Data_Zilionis2019_Lung.tar.gz?rlkey=kxqcfab6vaj043bvfkiepx2w9&dl=1",
        "metadata_url": "https://www.dropbox.com/scl/fi/0npq2qqtbap05mlzloqvd/Meta-data_Zilionis2019_Lung.tar.gz?rlkey=j9h9kck83ybltirr2t0cnudna&dl=1",
    },
}

# Bulk download options
BULK_DOWNLOADS = {
    "all_data": "https://www.dropbox.com/scl/fo/7wyp7x2irndqhn2drrc0c/AFAqeJ6K7_a6TqXb95UIc1E?rlkey=qjefasm7oha10f7yg74hpmz9q&dl=1",
    "all_metadata": "https://www.dropbox.com/scl/fo/8d0vi4pjnuf4gij3ll262/ADuPk5vj35q71RtELC8vz-c?rlkey=nz112f59ewp1v0rjd4mhjyd8r&dl=1",
}


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Download Curated Cancer Cell Atlas - Lung Cancer datasets",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("./LungAtlas"),
        help="Output directory for downloads (default: ./LungAtlas)",
    )
    parser.add_argument(
        "--studies",
        nargs="+",
        choices=list(DATASETS.keys()),
        help="Specific studies to download (e.g., Bischoff2021 Maynard2020)",
    )
    parser.add_argument(
        "--download-all",
        action="store_true",
        help="Download all available studies",
    )
    parser.add_argument(
        "--data-only",
        action="store_true",
        help="Download only data files (no metadata)",
    )
    parser.add_argument(
        "--metadata-only",
        action="store_true",
        help="Download only metadata files (no data)",
    )
    parser.add_argument(
        "--no-extract",
        action="store_true",
        help="Do not extract tar.gz files after download",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Skip files that already exist",
    )
    parser.add_argument(
        "--list-studies",
        action="store_true",
        help="List all available studies and exit",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging",
    )
    return parser.parse_args()


def list_studies() -> None:
    """Print all available studies."""
    print("\nAvailable Lung Cancer Studies in 3CA:")
    print("=" * 80)
    for study_id, info in DATASETS.items():
        print(f"\n{study_id}")
        print(f"  Full name: {info['full_name']}")
        print(f"  Description: {info['description']}")
    print("\n" + "=" * 80)
    print(f"Total: {len(DATASETS)} studies")
    print()


def format_bytes(size: int) -> str:
    """Format bytes to human-readable string."""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size < 1024.0:
            return f"{size:.1f} {unit}"
        size /= 1024.0
    return f"{size:.1f} PB"


def download_file(url: str, output_path: Path, resume: bool = False) -> bool:
    """Download a file from URL with progress tracking.

    Args:
        url: URL to download from
        output_path: Path to save the file
        resume: Skip if file already exists

    Returns:
        True if download succeeded, False otherwise
    """
    if resume and output_path.exists():
        logger.info(f"File already exists, skipping: {output_path.name}")
        return True

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = output_path.with_suffix(output_path.suffix + '.tmp')

    try:
        logger.info(f"Downloading {output_path.name}...")

        # Create request with user agent to avoid blocking
        req = Request(url, headers={'User-Agent': 'Mozilla/5.0'})

        with urlopen(req, timeout=30) as response:
            total_size = int(response.headers.get('content-length', 0))

            if total_size > 0:
                logger.info(f"File size: {format_bytes(total_size)}")

            downloaded = 0
            chunk_size = 8192
            start_time = time.time()
            last_log_time = start_time

            with open(temp_path, 'wb') as f:
                while True:
                    chunk = response.read(chunk_size)
                    if not chunk:
                        break

                    f.write(chunk)
                    downloaded += len(chunk)

                    # Log progress every 5 seconds
                    current_time = time.time()
                    if current_time - last_log_time >= 5.0:
                        if total_size > 0:
                            pct = (downloaded / total_size) * 100
                            elapsed = current_time - start_time
                            speed = downloaded / elapsed if elapsed > 0 else 0
                            logger.info(
                                f"Progress: {format_bytes(downloaded)} / {format_bytes(total_size)} "
                                f"({pct:.1f}%) - {format_bytes(speed)}/s"
                            )
                        else:
                            logger.info(f"Downloaded: {format_bytes(downloaded)}")
                        last_log_time = current_time

            # Move temp file to final location
            temp_path.rename(output_path)

            elapsed = time.time() - start_time
            avg_speed = downloaded / elapsed if elapsed > 0 else 0
            logger.info(
                f"Download complete: {output_path.name} "
                f"({format_bytes(downloaded)} in {elapsed:.1f}s, avg {format_bytes(avg_speed)}/s)"
            )
            return True

    except (URLError, HTTPError) as exc:
        logger.error(f"Download failed for {output_path.name}: {exc}")
        if temp_path.exists():
            temp_path.unlink()
        return False
    except Exception as exc:
        logger.exception(f"Unexpected error downloading {output_path.name}: {exc}")
        if temp_path.exists():
            temp_path.unlink()
        return False


def extract_tarball(archive_path: Path, extract_dir: Path) -> bool:
    """Extract a tar.gz file.

    Args:
        archive_path: Path to the tar.gz file
        extract_dir: Directory to extract to

    Returns:
        True if extraction succeeded, False otherwise
    """
    try:
        import tarfile

        logger.info(f"Extracting {archive_path.name}...")
        extract_dir.mkdir(parents=True, exist_ok=True)

        with tarfile.open(archive_path, 'r:gz') as tar:
            tar.extractall(path=extract_dir)

        logger.info(f"Extraction complete: {archive_path.name}")
        return True

    except Exception as exc:
        logger.exception(f"Failed to extract {archive_path.name}: {exc}")
        return False


def download_study(
    study_id: str,
    output_dir: Path,
    download_data: bool = True,
    download_metadata: bool = True,
    resume: bool = False,
    extract: bool = True,
) -> Dict[str, bool]:
    """Download data and/or metadata for a single study.

    Args:
        study_id: Study identifier
        output_dir: Base output directory
        download_data: Whether to download data files
        download_metadata: Whether to download metadata files
        resume: Skip files that already exist
        extract: Extract tar.gz files after download

    Returns:
        Dictionary with download results
    """
    if study_id not in DATASETS:
        logger.error(f"Unknown study: {study_id}")
        return {"data": False, "metadata": False}

    study_info = DATASETS[study_id]
    study_dir = output_dir / study_id
    results = {"data": None, "metadata": None}

    logger.info(f"\n{'=' * 80}")
    logger.info(f"Downloading {study_info['full_name']}")
    logger.info(f"Description: {study_info['description']}")
    logger.info(f"Output: {study_dir}")
    logger.info(f"{'=' * 80}\n")

    # Download data
    if download_data:
        data_file = study_dir / f"Data_{study_id}_Lung.tar.gz"
        results["data"] = download_file(study_info["data_url"], data_file, resume)

        if results["data"] and extract:
            extract_dir = study_dir / "data"
            extract_tarball(data_file, extract_dir)

    # Download metadata
    if download_metadata:
        metadata_file = study_dir / f"Meta-data_{study_id}_Lung.tar.gz"
        results["metadata"] = download_file(study_info["metadata_url"], metadata_file, resume)

        if results["metadata"] and extract:
            extract_dir = study_dir / "metadata"
            extract_tarball(metadata_file, extract_dir)

    return results


def main() -> None:
    """Main entry point."""
    args = parse_args()

    # Setup logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # List studies and exit
    if args.list_studies:
        list_studies()
        sys.exit(0)

    # Validate arguments
    if not args.download_all and not args.studies:
        logger.error("Must specify either --download-all or --studies")
        logger.info("Use --list-studies to see available studies")
        sys.exit(1)

    if args.data_only and args.metadata_only:
        logger.error("Cannot specify both --data-only and --metadata-only")
        sys.exit(1)

    # Determine what to download
    download_data = not args.metadata_only
    download_metadata = not args.data_only

    # Determine which studies to download
    if args.download_all:
        studies_to_download = list(DATASETS.keys())
    else:
        studies_to_download = args.studies

    # Print summary
    print("\n" + "=" * 80)
    print("Download Summary")
    print("=" * 80)
    print(f"Output directory: {args.output_dir}")
    print(f"Studies to download: {len(studies_to_download)}")
    print(f"Download data: {download_data}")
    print(f"Download metadata: {download_metadata}")
    print(f"Extract after download: {not args.no_extract}")
    print(f"Resume mode: {args.resume}")
    print("=" * 80 + "\n")

    # Download studies
    start_time = time.time()
    success_count = 0
    failure_count = 0

    for study_id in studies_to_download:
        results = download_study(
            study_id,
            args.output_dir,
            download_data=download_data,
            download_metadata=download_metadata,
            resume=args.resume,
            extract=not args.no_extract,
        )

        # Count successes
        if download_data and results.get("data"):
            success_count += 1
        elif download_data and results.get("data") is False:
            failure_count += 1

        if download_metadata and results.get("metadata"):
            success_count += 1
        elif download_metadata and results.get("metadata") is False:
            failure_count += 1

    elapsed = time.time() - start_time

    # Print final summary
    print("\n" + "=" * 80)
    print("Download Complete")
    print("=" * 80)
    print(f"Total time: {elapsed / 60:.1f} minutes")
    print(f"Studies processed: {len(studies_to_download)}")
    print(f"Successful downloads: {success_count}")
    print(f"Failed downloads: {failure_count}")
    print(f"Output directory: {args.output_dir}")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    main()
