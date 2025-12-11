import os
import glob
import requests
import re
import sys

# Configuration
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
POSTS_DIR = os.path.join(BASE_DIR, "../content/posts")
DATA_DIR = os.path.join(BASE_DIR, "../data")
KEY_FILE = os.path.join(DATA_DIR, "UNSPLASH_ACCESS_KEY")

def get_access_key():
    try:
        with open(KEY_FILE, 'r') as f:
            return f.read().strip()
    except FileNotFoundError:
        print(f"Warning: Unsplash access key file not found at {KEY_FILE}")
        return None

UNSPLASH_ACCESS_KEY = get_access_key()

def get_random_photo():
    if not UNSPLASH_ACCESS_KEY:
        print("Error: UNSPLASH_ACCESS_KEY not found. Please ensure data/UNSPLASH_ACCESS_KEY exists.")
        return None
    
    url = "https://api.unsplash.com/photos/random?orientation=landscape"
    headers = {"Authorization": f"Client-ID {UNSPLASH_ACCESS_KEY}"}
    
    try:
        response = requests.get(url, headers=headers)
        response.raise_for_status()
        data = response.json()
        return data["urls"]["regular"]
    except Exception as e:
        print(f"Failed to fetch image from Unsplash: {e}")
        return None

def process_file(filepath, force=False):
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Check if front matter exists
    if not content.startswith("+++"):
        return

    # Extract front matter
    parts = content.split("+++", 2)
    if len(parts) < 3:
        return
    
    front_matter = parts[1]
    body = parts[2]
    
    # Check if feature_image is already set to a specific static image
    # We want to replace generic random URLs
    match = re.search(r'feature_image\s*=\s*"(.*?)"', front_matter)
    if match:
        current_image = match.group(1)
        # If it's a specific Unsplash photo ID (usually long) or a local path, keep it.
        # If it's a generic randomizer, replace it.
        is_generic = "random" in current_image or "source.unsplash.com" in current_image or "dmoe.cc" in current_image
        
        if not is_generic and not force:
            print(f"Skipping {filepath}: feature_image already set to {current_image}")
            return
        
        if is_generic or force:
            if force:
                print(f"Force overwriting image in {filepath}...")
            else:
                print(f"Overwriting generic image in {filepath}...")
            # Remove the existing line so we can add the new one
            front_matter = re.sub(r'feature_image\s*=\s*".*?"\n?', '', front_matter)

    print(f"Processing {filepath}...")
    image_url = get_random_photo()
    
    if image_url:
        # Add feature_image to [extra] section or create it
        if "[extra]" in front_matter:
            front_matter = front_matter.replace("[extra]", f"[extra]\nfeature_image = \"{image_url}\"")
        else:
            front_matter += f"\n[extra]\nfeature_image = \"{image_url}\"\n"
        
        new_content = f"+++{front_matter}+++{body}"
        
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(new_content)
        print(f"Updated {filepath} with image: {image_url}")
    else:
        print(f"Skipping {filepath}: Could not get image.")

def main():
    if not UNSPLASH_ACCESS_KEY:
        print("Please ensure the Unsplash Access Key is in data/UNSPLASH_ACCESS_KEY.")
        return

    force = "--force" in sys.argv

    files = glob.glob(os.path.join(POSTS_DIR, "*.md"))
    for filepath in files:
        if "_index.md" in filepath:
            continue
        process_file(filepath, force)

if __name__ == "__main__":
    main()
