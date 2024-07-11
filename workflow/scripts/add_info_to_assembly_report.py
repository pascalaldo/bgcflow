import sys
import json

def add_info(json_path):
    try:
        with open(json_path, "r") as f:
            d = json.load(f)
        genus, species = d["organism"]["organismName"].split(" ", 1)
        strain = d["organism"].get("infraspecificNames", {}).get("strain", "")
        taxid = d["organism"]["taxId"]
        d["genus"] = genus
        d["species"] = species
        d["strain"] = strain
        d["taxid"] = taxid
        with open(json_path, "w") as f:
            json.dump(d, f)
    except:
        print(f"Warning: Could not load '{json_path}', not adding info.")

if __name__ == "__main__":
    add_info(sys.argv[1])