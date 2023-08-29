# main script
from SHPB_GUI import SHPB_GUI
import tkinter as tk

if __name__ == '__main__':
    #run the GUI
    root = tk.Tk()
    gui = SHPB_GUI(root)
    root.mainloop()

#ideas: 
#expand every csv in the GUI and report the quality of dynamic equilibrium (with warning for bad results), add button for each that shows the graph
#option to manually alter equilibrium (similar to matlab script) 
#customizeable preset buttons for aluminum and steel bars
#least squares inspired equilibrium verification (maybe faster?)
#integrate and subtract equilibrium verification (slower but more accurate?)
#option to upload separate csvs for incident and transmission bar data (convention with other software)
#parallelization to do several csvs at once (speed increase)
#look into generalizing for tension, torsion
#look into overlap data (ask andrew for sample overlaps?)