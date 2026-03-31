import os
# Must set this BEFORE importing marker modules
# os.environ["TORCH_DEVICE"] = "cpu"
import re
from marker.converters.pdf import PdfConverter
from marker.models import create_model_dict
# The new structure for 2026
def sanitize_filename(name):
    return re.sub(r'[\\/*?:"<>|]', "", name).strip()

def convert_folder_to_markdown(input_folder, output_folder="./Converted_Markdown"):
    print("--- Initializing AI Models (CPU Mode) ---")
    # create_model_dict() now handles the device based on the environment variable
    converter = PdfConverter(
        artifact_dict=create_model_dict(),
    )
    
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    pdf_files = [f for f in os.listdir(input_folder) if f.lower().endswith('.pdf')]
    
    if not pdf_files:
        print(f"No PDFs found in {input_folder}")
        return

    print(f"Converting {len(pdf_files)} papers on CPU. This may take a while...\n")

    for i, pdf_name in enumerate(pdf_files, 1):
        input_path = os.path.join(input_folder, pdf_name)
        base_name = sanitize_filename(os.path.splitext(pdf_name)[0])
        
        # Create a subfolder for each paper's assets (images + md)
        paper_dir = os.path.join(output_folder, base_name)
        if not os.path.exists(paper_dir):
            os.makedirs(paper_dir)
            
        md_file_path = os.path.join(paper_dir, f"{base_name}.md")
        
        if os.path.exists(md_file_path):
            print(f"[{i}/{len(pdf_files)}] Skipping: {pdf_name}")
            continue

        print(f"[{i}/{len(pdf_files)}] Converting: {pdf_name}...")
        
        try:
            # Execute conversion
            rendered = converter(input_path)
            
            # 2. SAVE MANUALLY (Fixes the 'save_all' error)
            # Access the markdown text via the .markdown attribute
            with open(md_file_path, "w", encoding="utf-8") as f:
                f.write(rendered.markdown)
            
            # Save images if they exist
            if rendered.images:
                image_dir = os.path.join(paper_dir, "images")
                if not os.path.exists(image_dir):
                    os.makedirs(image_dir)
                for img_name, img_data in rendered.images.items():
                    with open(os.path.join(image_dir, img_name), "wb") as img_f:
                        img_f.write(img_data)
            
            print(f"   Successfully saved to {paper_dir}")
            
        except Exception as e:
            print(f"Error on {pdf_name}: {e}")

    print(f"\n--- Batch processing complete! ---")

if __name__ == "__main__":
    # Ensure this path exists
    target_folder = "./Tepper_Gill_Papers"
    if os.path.exists(target_folder):
        convert_folder_to_markdown(target_folder)
    else:
        print(f"Target folder '{target_folder}' not found.")