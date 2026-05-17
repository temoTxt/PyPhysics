"""PDF → Markdown converter for the PyPhysics paper corpus.

Wraps `marker-pdf` to convert every PDF in an input directory into a
per-paper subdirectory of markdown + extracted images. Used both for the
original Gill corpus (`Tepper_Gill_Papers/` → `Converted_Markdown/`) and
for the History campaign's historical papers (`Historical_Papers/Primary/`
→ `Historical_Converted_Markdown/Primary/`, and same for `Retrospective/`).

Usage:
    uv run python Roadmapping/parse_papers.py                                 # default Gill corpus
    uv run python Roadmapping/parse_papers.py --input <in> --output <out>     # arbitrary tree
    uv run python Roadmapping/parse_papers.py --input <in> --output <out> --no-skip-existing

Run from the repository root so the default paths resolve correctly. The
conversion runs on CPU by default; export `TORCH_DEVICE=cpu` if torch is
trying to grab a GPU you don't want it to use.
"""

import argparse
import os
import re
import sys

DEFAULT_INPUT = "./Roadmapping/Tepper_Gill_Papers"
DEFAULT_OUTPUT = "./Roadmapping/Converted_Markdown"


def sanitize_filename(name: str) -> str:
    return re.sub(r'[\\/*?:"<>|]', "", name).strip()


def convert_folder_to_markdown(input_folder: str, output_folder: str, skip_existing: bool = True) -> int:
    """Convert every PDF in input_folder to markdown under output_folder/<paper>/<paper>.md.

    Returns the number of papers successfully converted (skips not counted).
    """
    # Imported here so --help works without the heavy ML deps installed.
    from marker.converters.pdf import PdfConverter
    from marker.models import create_model_dict

    print("--- Initializing AI Models (CPU Mode) ---")
    converter = PdfConverter(artifact_dict=create_model_dict())

    os.makedirs(output_folder, exist_ok=True)

    pdf_files = [f for f in os.listdir(input_folder) if f.lower().endswith(".pdf")]
    if not pdf_files:
        print(f"No PDFs found in {input_folder}")
        return 0

    print(f"Converting {len(pdf_files)} papers on CPU. This may take a while...\n")
    converted = 0
    for i, pdf_name in enumerate(pdf_files, 1):
        input_path = os.path.join(input_folder, pdf_name)
        base_name = sanitize_filename(os.path.splitext(pdf_name)[0])
        paper_dir = os.path.join(output_folder, base_name)
        os.makedirs(paper_dir, exist_ok=True)
        md_file_path = os.path.join(paper_dir, f"{base_name}.md")

        if skip_existing and os.path.exists(md_file_path):
            print(f"[{i}/{len(pdf_files)}] Skipping (already converted): {pdf_name}")
            continue

        print(f"[{i}/{len(pdf_files)}] Converting: {pdf_name}...")
        try:
            rendered = converter(input_path)
            with open(md_file_path, "w", encoding="utf-8") as f:
                f.write(rendered.markdown)

            if rendered.images:
                image_dir = os.path.join(paper_dir, "images")
                os.makedirs(image_dir, exist_ok=True)
                for img_name, img_data in rendered.images.items():
                    with open(os.path.join(image_dir, img_name), "wb") as img_f:
                        img_f.write(img_data)

            print(f"   Successfully saved to {paper_dir}")
            converted += 1
        except Exception as e:
            print(f"Error on {pdf_name}: {e}")

    print(f"\n--- Batch processing complete: {converted} converted, "
          f"{len(pdf_files) - converted} skipped or failed. ---")
    return converted


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", default=DEFAULT_INPUT,
                        help=f"directory of PDFs to convert (default: {DEFAULT_INPUT})")
    parser.add_argument("--output", default=DEFAULT_OUTPUT,
                        help=f"output directory for per-paper subdirs (default: {DEFAULT_OUTPUT})")
    parser.add_argument("--skip-existing", dest="skip_existing", action="store_true", default=True,
                        help="skip PDFs whose markdown already exists (default)")
    parser.add_argument("--no-skip-existing", dest="skip_existing", action="store_false",
                        help="re-convert every PDF even if markdown already exists")
    args = parser.parse_args()

    if not os.path.isdir(args.input):
        print(f"Input directory not found: {args.input}", file=sys.stderr)
        return 1

    convert_folder_to_markdown(args.input, args.output, skip_existing=args.skip_existing)
    return 0


if __name__ == "__main__":
    sys.exit(main())
