import tkinter as tk
from tkinter import filedialog, messagebox
from padelpy import padeldescriptor
import os

def browse_file():
    """Open file browser and select .sdf file."""
    filename = filedialog.askopenfilename(
        title="Select SDF File",
        filetypes=[("SDF files", "*.sdf"), ("All files", "*.*")]
    )
    if filename:
        sdf_entry.delete(0, tk.END)
        sdf_entry.insert(0, filename)

def run_job():
    """Run PaDEL descriptor calculation based on user choices."""
    sdf_file = sdf_entry.get()
    if not sdf_file or not os.path.exists(sdf_file):
        messagebox.showerror("Error", "Please select a valid SDF file.")
        return
    
    output_csv = os.path.splitext(sdf_file)[0] + "_descriptors.csv"
    
    # Set user options
    use_2d = desc_var.get()
    use_3d = desc3d_var.get()
    use_fp = fp_var.get()
    
    
    if use_2d:
       d_2d=True
    else:
       d_2d=False
    if use_3d:
       d_3d=True
    else:
       d_3d=False
    if use_fp:
       use_fp=True
    else:
       use_fp=False


       
    # Build arguments for PaDEL
    padeldescriptor(
        mol_dir=sdf_file,
        d_file=output_csv,          # output CSV file
        d_2d=d_2d,                  # calculate 2D descriptors
        d_3d=d_3d,               # calculate 3D descriptors if coords exist
        fingerprints=use_fp,  # optional: calculate fingerprints
        threads=4,
        removesalt=True,
        standardizenitro=True,
        retainorder=True,           # keep order same as input
        headless=True,
        log=True
    )

    msg = f"Descriptors saved to:\n{output_csv}\n\n"
    msg += "Options selected:\n"
    if use_2d: msg += "- 1D/2D descriptors\n"
    if use_3d: msg += "- 3D descriptors (⚠ PaDEL does not support true 3D!)\n"
    if use_fp: msg += "- Fingerprints\n"

    messagebox.showinfo("Job Completed", msg)

# ---------------- TKINTER GUI ---------------- #
root = tk.Tk()
root.title("PaDEL Descriptor Calculator")

# File browser
tk.Label(root, text="SDF File:").grid(row=0, column=0, padx=5, pady=5, sticky="e")
sdf_entry = tk.Entry(root, width=40)
sdf_entry.grid(row=0, column=1, padx=5, pady=5)
browse_btn = tk.Button(root, text="Browse", command=browse_file)
browse_btn.grid(row=0, column=2, padx=5, pady=5)

# Options
desc_var = tk.BooleanVar(value=True)
desc3d_var = tk.BooleanVar(value=False)
fp_var = tk.BooleanVar(value=True)

tk.Checkbutton(root, text="1D/2D Descriptors", variable=desc_var).grid(row=1, column=0, columnspan=2, sticky="w", padx=10)
tk.Checkbutton(root, text="3D Descriptors", variable=desc3d_var).grid(row=2, column=0, columnspan=2, sticky="w", padx=10)
tk.Checkbutton(root, text="Fingerprints", variable=fp_var).grid(row=3, column=0, columnspan=2, sticky="w", padx=10)

# Run Job button
run_btn = tk.Button(root, text="Run Job", command=run_job, bg="green", fg="white")
run_btn.grid(row=4, column=1, pady=10)

root.mainloop()
