import numpy as np
import matplotlib.pyplot as plt
import os
import tkinter as tk
from tkinter import filedialog
from tkinter import ttk
from SHPB_analysis import SHPB_analysis

class SHPB_GUI:
    def __init__(self, master):
        self.master = master
        self.master.title("Split-Hopkinson Pressure Bar Analysis")
        self.master.geometry("350x450")

        #bar and sample properties
        self.young_modulus = tk.DoubleVar(value=210.0e+09)
        self.density = tk.DoubleVar(value=8100)
        self.rate = tk.DoubleVar(value=25000000)
        self.sample_area = tk.DoubleVar(value=0.01)
        self.bar_area = tk.DoubleVar()
        self.bar_diameter = tk.DoubleVar(value=0.01905)
        self.voltage_to_strain_factor = tk.DoubleVar(value=0.001)
        self.sample_mass = tk.DoubleVar(value=0.000079)
        self.sample_thickness =  tk.DoubleVar(value=0.00054)
        self.folder = tk.StringVar()
        
        self.create_widgets()

    def create_widgets(self):
        
        #GUI setup
        mainframe = ttk.Frame(self.master, padding="3 3 12 12")
        mainframe.grid(column=0, row=0, sticky=(tk.N, tk.W, tk.E, tk.S))
        mainframe.columnconfigure(0, weight=1)
        mainframe.rowconfigure(0, weight=1)

        #input boxes for bar and sample properties
        self.young_modulus_entry = ttk.Entry(mainframe, textvariable=self.young_modulus)
        self.young_modulus_entry.grid(column=2, row=1, sticky=(tk.W, tk.E))

        self.density_entry = ttk.Entry(mainframe, textvariable=self.density)
        self.density_entry.grid(column=2, row=2, sticky=(tk.W, tk.E))

        self.rate_entry = ttk.Entry(mainframe, textvariable=self.rate)
        self.rate_entry.grid(column=2, row=3, sticky=(tk.W, tk.E))

        self.bar_diameter_entry = ttk.Entry(mainframe, textvariable=self.bar_diameter)
        self.bar_diameter_entry.grid(column=2, row=4, sticky=(tk.W, tk.E))

        self.voltage_to_strain_factor_entry = ttk.Entry(mainframe, textvariable=self.voltage_to_strain_factor)
        self.voltage_to_strain_factor_entry.grid(column=2, row=5, sticky=(tk.W, tk.E))

        self.sample_mass_entry = ttk.Entry(mainframe, textvariable=self.sample_mass)
        self.sample_mass_entry.grid(column=2, row=6, sticky=(tk.W, tk.E))

        self.sample_thickness_entry = ttk.Entry(mainframe, textvariable=self.sample_thickness)
        self.sample_thickness_entry.grid(column=2, row=7, sticky=(tk.W, tk.E))

        #input button to allow folder selection
        self.folder_button = ttk.Button(mainframe, text="Select Folder", command=lambda: self.folder.set(filedialog.askdirectory()))
        self.folder_button.grid(column=2, row=8, sticky=(tk.W, tk.E))

        #run button
        self.run_button = ttk.Button(mainframe, text="Run", command=self.run_analysis_callback)
        self.run_button.grid(column=2, row=9, sticky=(tk.W, tk.E))

        #labels
        ttk.Label(mainframe, text="Bar Young's Modulus (Pa)").grid(column=1, row=1, sticky=tk.E)
        ttk.Label(mainframe, text="Bar Density (kg/m^3)").grid(column=1, row=2, sticky=tk.E)
        ttk.Label(mainframe, text="Sampling Rate (Hz)").grid(column=1, row=3, sticky=tk.E)
        ttk.Label(mainframe, text="Bar Diameter (m)").grid(column=1, row=4, sticky=tk.E)
        ttk.Label(mainframe, text="Voltage to Strain Factor (V to ue)").grid(column=1, row=5, sticky=tk.E)
        ttk.Label(mainframe, text="Sample Mass (kg)").grid(column=1, row=6, sticky=tk.E)
        ttk.Label(mainframe, text="Sample Thickness (m)").grid(column=1, row=7, sticky=tk.E)

        for child in mainframe.winfo_children(): 
            child.grid_configure(padx=5, pady=5)

        self.young_modulus_entry.focus()
        self.master.bind('<Return>', self.run_analysis_callback)

        #store the selected sample shape
        self.sample_shape = tk.StringVar(value="circular")

        #drop-down menu for selecting the sample shape
        self.sample_shape_option_menu = ttk.OptionMenu(mainframe, self.sample_shape, "rectangular", "rectangular", "circular")
        self.sample_shape_option_menu.grid(column=2, row=10, sticky=(tk.W, tk.E))
        ttk.Label(mainframe, text="Sample Shape").grid(column=1, row=10, sticky=tk.E)

        #sample length, width, and diameter
        self.sample_length = tk.StringVar(value=0.01103)
        self.sample_width = tk.StringVar(value=0.01163)
        self.sample_diameter = tk.StringVar(value=0.01)

        #input boxes for the sample length and width
        self.sample_length_entry = ttk.Entry(mainframe, textvariable=self.sample_length)
        self.sample_length_entry.grid(column=2, row=11, sticky=(tk.W, tk.E))
        ttk.Label(mainframe, text="Sample Length (m)").grid(column=1, row=11, sticky=tk.E)

        self.sample_width_entry = ttk.Entry(mainframe, textvariable=self.sample_width)
        self.sample_width_entry.grid(column=2, row=12, sticky=(tk.W, tk.E))
        ttk.Label(mainframe, text="Sample Width (m)").grid(column=1, row=12, sticky=tk.E)

        #input box for the sample diameter
        self.sample_diameter_entry = ttk.Entry(mainframe, textvariable=self.sample_diameter)
        self.sample_diameter_entry.grid(column=2, row=13, sticky=(tk.W, tk.E))
        ttk.Label(mainframe, text="Sample Diameter (m)").grid(column=1, row=13, sticky=tk.E)

        #create a function that will hide or show the widgets based on the selected option
        self.sample_shape.trace('w', self.on_sample_shape_change)
        self.on_sample_shape_change()
        
    def on_sample_shape_change(self, *args):
        if self.sample_shape.get() == 'rectangular':
            self.sample_length_entry.grid(column=2, row=11, sticky=(tk.W, tk.E))
            self.sample_width_entry.grid(column=2, row=12, sticky=(tk.W, tk.E))
            self.sample_diameter_entry.grid_remove()
        else:
            self.sample_length_entry.grid_remove()
            self.sample_width_entry.grid_remove()
            self.sample_diameter_entry.grid(column=2, row=13, sticky=(tk.W, tk.E))
    
    def run_analysis_callback(self):
        print("running")

        if self.sample_shape.get() == "circular":
            try:
                self.sample_area = np.pi*(float(self.sample_diameter.get())**2)/4
            except ValueError:
                tk.messagebox.showerror("Error", "Invalid value entered for sample diameter. Please enter a numeric value.")
                return
            
        elif self.sample_shape.get() == "rectangular":
            try:
                self.sample_area = float(self.sample_length.get()) * float(self.sample_width.get())
            except ValueError:
                tk.messagebox.showerror("Error", "Invalid value entered for sample length or width. Please enter a numeric value.")
                return

        self.bar_area = np.pi*(float(self.bar_diameter.get())**2)/4
    
        csv_files = [f for f in os.listdir(self.folder.get()) if f.endswith('.csv') or f.endswith('.CSV')]

        for file in csv_files:
            file_path = os.path.join(self.folder.get(), file)

            SHPB_analysis.run_analysis(
                self, 
                self.young_modulus.get(), 
                self.density.get(), 
                self.rate.get(), 
                self.sample_area, 
                self.bar_area, 
                self.voltage_to_strain_factor.get(), 
                self.sample_mass.get(), 
                self.sample_thickness.get(),
                file_path
            )
        