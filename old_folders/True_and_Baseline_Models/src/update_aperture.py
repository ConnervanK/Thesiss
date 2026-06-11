 import sys
import old_folders.True_and_Baseline_Models.src.helper as helper

def main():
    if len(sys.argv) < 3:
        print("Usage: python update_aperture.py <filename> <new_aperture_in_meters> [fixed_reference] [material_name]")
        print("fixed_reference options: top (default), bottom, center")
        print("material_name: default 'water'")
        print("Example: python update_aperture.py process_model.in 0.005 top water")
        return

    filename = sys.argv[1]
    
    try:
        new_aperture = float(sys.argv[2])
    except ValueError:
        print("Error: new_aperture must be a number.")
        return

    fixed_reference = 'top'
    if len(sys.argv) > 3:
        fixed_reference = sys.argv[3]
    
    material_name = 'water'
    if len(sys.argv) > 4:
        material_name = sys.argv[4]

    print(f"Updating {filename} with aperture {new_aperture} m.")
    print(f"  Fixed reference: {fixed_reference}")
    print(f"  Material: {material_name}")
    
    # We pass material_name now
    helper.adjust_fracture_aperture(filename, new_aperture, fixed_reference, material_name)

if __name__ == "__main__":
    main()
