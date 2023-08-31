import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.widgets import Button
from scipy.integrate import cumtrapz
from scipy.signal import butter, filtfilt
from scipy.optimize import minimize
from scipy.ndimage.interpolation import shift



class SHPB_analysis:

    def run_analysis(self, 
        young_modulus, 
        density, 
        rate, 
        sample_area, 
        bar_area, 
        voltage_to_strain_factor, 
        sample_mass, 
        sample_thickness,
        file_path
    ):

        data = pd.read_csv(file_path, header=None)

        #trim and cast the data
        data = data.iloc[3:,1:] #trim header and first 3 rows
        data = data.astype(float)

        #incident, transmission, and time dataframes
        data["incident_bar"] = data.iloc[:,0]
        data["transmitted_bar"] = data.iloc[:,1]
        data['time'] = np.arange(len(data)) * (1/rate)  # add time column based on sampling rate

        if np.max(data["incident_bar"]) > 5: #check scale
            data["incident_bar"] /= 1000

        if np.max(data["transmitted_bar"]) > 5: #check scale, make better solution later
            data["transmitted_bar"] /= 1000

        #calculate initial bias
        incident_bias = np.sum(data["incident_bar"].head(int(np.round(0.02*len(data["incident_bar"]))))/np.round(0.02*len(data["incident_bar"])))
      
        transmission_bias = np.sum(data["transmitted_bar"].head(int(np.round(0.02*len(data["transmitted_bar"]))))/np.round(0.02*len(data["transmitted_bar"])))
     

        #correct for bias
        data["incident_bar"] -= incident_bias
        data["transmitted_bar"] -= transmission_bias

        #plot signal
        plt.plot(data["time"].values*1e6, data["transmitted_bar"].values, label='Transmitted Bar')
        plt.plot(data["time"].values*1e6, data["incident_bar"].values, label='Incident Bar')
        plt.xlabel('Time (us)')
        plt.ylabel('Amplitude (V)')
        plt.title('SHPB Analysis')
        plt.legend()
        plt.show()

        #find incident wave 

        global reflected_wave, transmitted_wave, incident_wave

        #NORMALIZE THIS AND GO BY A PERCENTAGE
    
        #define indices 
        start_index_inc = None
        end_index_inc = None
        
        #normalize waves
        trans_min = min(data["transmitted_bar"].values) #when adapting for tension, use max or use more robust wave identification method
        inc_min = min(data["incident_bar"].values)

        norm_transmitted_bar = -1*data["transmitted_bar"].values/trans_min
        norm_incident_bar = -1*data["incident_bar"]/inc_min

        #filter waves

        # Define the cutoff frequency and order of the filter
        cutoff = 0.1
        order = 6

        # Create the Butterworth filter coefficients
        b, a = butter(order, cutoff, btype='low', analog=False, output='ba')

        # Apply the filter to the signal
        norm_transmitted_bar_filtered = filtfilt(b, a, norm_transmitted_bar)
        norm_incident_bar_filtered = filtfilt(b, a, norm_incident_bar)

        # Plot
        plt.plot(data["time"].values, norm_transmitted_bar_filtered, label="Filtered Transmitted")
        plt.plot(data["time"].values, norm_incident_bar_filtered, label="Filtered Incident")
        plt.xlabel('Time (s)')
        plt.ylabel('Amplitude')
        plt.legend()
        plt.show()

        for i, value in enumerate(norm_incident_bar_filtered):
            if value < -0.1: #finding the wave
                if start_index_inc is None:
                    start_index_inc = i      
            if value > -0.1:
                if start_index_inc is not None:
                    end_index_inc = i
                    if abs(start_index_inc - end_index_inc) < 20:
                        start_index_inc = None
                    else:
                        break
        #print(start_index)
        #print(end_index)

        for i in range(start_index_inc, 0, -1): # iterate backwards starting from start_index
            if norm_incident_bar_filtered[i] > 0.0: #can get away with 0.0 because of bias correction
                start_index_inc = i
                break
        #print(start_index)

        for i in range(end_index_inc, len(norm_incident_bar_filtered)): # iterate forwards starting from end_index
            if norm_incident_bar_filtered[i] > -0.02: #ask damian for look at brittle samples to see how much i can get away with here
                end_index_inc = i + round(0.02*i)
                break
        #print(end_index)

        #create incident wave array
        incident_wave = data["incident_bar"].values[start_index_inc:end_index_inc].astype(dtype='double')

        #finding transmission wave
        start_index_trans = None
        end_index_trans = None

        for i, value in enumerate(norm_transmitted_bar_filtered):
            if value < -0.1: #finding the wave
                if start_index_trans is None:
                    start_index_trans = i      
            if value > -0.1:
                if start_index_trans is not None:
                    end_index_trans = i
                    if abs(start_index_trans - end_index_trans) < 20:
                        start_index_trans = None
                    else:
                        break
        #print(start_index_trans)
        #print(end_index_trans)

        for i in range(start_index_trans, 0, -1): # iterate backwards starting from start_index
            if norm_transmitted_bar_filtered[i] > -0.005: 
                start_index_trans = i - round(0.01*i)
                break
        #print(start_index_trans)

        end_index_trans = start_index_trans + len(incident_wave)

        #create transmitted wave array

        transmitted_wave = data["transmitted_bar"].values[start_index_trans:end_index_trans].astype(dtype='double')

        #finding reflected wave
        start_index_ref = None
        end_index_ref = None

        for i, value in enumerate(norm_incident_bar_filtered):
            if value > 0.1: #finding the wave
                if start_index_ref is None:
                    start_index_ref = i      
            if value < 0.1:
                if start_index_ref is not None:
                    end_index_ref = i
                    if abs(start_index_ref - end_index_ref) < 20:
                        start_index_ref = None
                    else:
                        break
        #print(start_index_ref)
        #print(end_index_ref)

        for i in range(start_index_ref, 0, -1): # iterate backwards starting from start_index
            #print(i)
            if norm_incident_bar_filtered[i] < 0.005:
                start_index_ref = i - round(0.01*i)
                break
        #print(start_index_ref)

        end_index_ref = start_index_ref + len(incident_wave)

        #create reflected wave array
        reflected_wave = data["incident_bar"].values[start_index_ref:end_index_ref].astype(dtype='double')

        #convert from voltage to strain (sorry for abusive notation)
        incident_wave *= voltage_to_strain_factor
        transmitted_wave *= voltage_to_strain_factor
        reflected_wave *= voltage_to_strain_factor

        time_window_end = end_index_inc - start_index_inc

        # plot the incident wave data
        plt.plot(data['time'].values[:time_window_end]*1e6, incident_wave, label='Incident Wave')
        plt.plot(data['time'].values[:time_window_end]*1e6, transmitted_wave, label='Transmitted Wave')
        plt.plot(data['time'].values[:time_window_end]*1e6, reflected_wave, label='Reflected Wave')
        plt.xlabel('Time (us)')
        plt.ylabel('Amplitude')
        plt.title('SHPB Analysis')
        plt.legend()
        plt.show()

        #####################################
        #dynamic stress equilibrium
        #####################################

        #cost function for force balance

        print(len(incident_wave))

        def cost_function(x, incident_wave, reflected_wave, transmitted_wave, shift_penalty=1e-8):
            shift1, shift2 = x
            reflected_wave_shifted = shift(reflected_wave, shift1, order=3, mode="constant", cval=0)
            transmitted_wave_shifted = shift(transmitted_wave, shift2, order=3, mode="constant", cval=0)

            # trans = ref + inc
            difference = incident_wave + reflected_wave_shifted - transmitted_wave_shifted

            # calculate the sum of squared differences
            sum_of_squared_differences = np.sum(difference**2)

            # calculate penalty for large shifts; normalize by power of 2/3 of wave size (increases a lot to start, increased wave size decreases importance of this)
            shift_penalty = shift_penalty/(len(incident_wave)**(2/3)) * (shift1**2 + shift2**2)

            # combine the cost and penalty terms
            total_cost = sum_of_squared_differences + shift_penalty

            return total_cost
        

        initial_guess = [0, 0]
        result = minimize(cost_function, initial_guess, args=(incident_wave, reflected_wave, transmitted_wave), method="powell") #cobyla or powell are best
        optimal_shift1, optimal_shift2 = result.x
        
        #initial_guess = [optimal_shift1, optimal_shift2]
        #result = minimize(cost_function, initial_guess, args=(incident_wave, reflected_wave, transmitted_wave), method="cobyla")
        #optimal_shift1, optimal_shift2 = result.x

        #from scipy.optimize import direct
        #bounds = [(-5000,5000),(-5000,5000)]  # define bounds for the shift parameters
        #result = direct(cost_function, bounds=bounds, args=(incident_wave, reflected_wave, transmitted_wave))
        #optimal_shift1, optimal_shift2 = result.x

        #initialize shift values
        global shift1, shift2
        shift1 = optimal_shift1
        shift2 = optimal_shift2

        reflected_wave = shift(reflected_wave, shift1, order=3, mode="nearest", cval=0)
        transmitted_wave = shift(transmitted_wave, shift2, order=3, mode="nearest", cval=0)
        incident_end = -1*young_modulus*bar_area*(incident_wave + reflected_wave)
        transmission_end = -1*young_modulus*bar_area*(transmitted_wave)
        average_force = (incident_end+transmission_end)/2
        

        #create figure and axes for both plots
        fig, axs = plt.subplots(1, 2, figsize=(10, 5))
        plt.subplots_adjust(left=0.25, bottom=0.25)

        #create first plot of force
        axs[0].plot(data['time'].values[:time_window_end]*1e6, incident_end, label='Incident End', color = 'b')
        axs[0].plot(data['time'].values[:time_window_end]*1e6, transmission_end, label='Transmission End', color = 'r')
        axs[0].plot(data['time'].values[:time_window_end]*1e6, average_force, label='Average', color = 'g')
        axs[0].set_xlabel('Time (us)')
        axs[0].set_ylabel('Force (N)')
        axs[0].set_title('SHPB Analysis')
        axs[0].legend()

        #create second plot of waves
        peak_index = reflected_wave.argmax()
        peak_time = data['time'].values[peak_index] * 1e6
        axs[1].plot(data['time'].values[:time_window_end]*1e6, incident_wave, label='Incident Wave', color = 'b')
        axs[1].plot(data['time'].values[:time_window_end]*1e6, reflected_wave, label='Reflected Wave', color = 'g')
        axs[1].axvline(peak_time, linestyle='--', color='grey')
        axs[1].plot(data['time'].values[:time_window_end]*1e6, transmitted_wave, label='Transmitted Wave', color = 'r')
        axs[1].set_xlabel('Time (us)')
        #axs[1].set_ylabel('Amplitude')
        axs[1].set_title('Waves')
        axs[1].legend()

        #create sliders for adjusting shifts
        axcolor = 'lightgoldenrodyellow'
        ax_shift1 = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
        ax_shift2 = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
        shift1_slider = Slider(ax_shift1, 'Reflected Shift', -500,500, valinit=0, valstep=1)
        shift2_slider = Slider(ax_shift2, 'Transmitted Shift', -500,500, valinit=0, valstep=1)

        #function to update plots when sliders are adjusted
        def update(val):
            global shift1, shift2
            global reflected_wave, transmitted_wave, incident_wave
            shift1 = shift1_slider.val
            shift2 = shift2_slider.val
            reflected_wave_new = shift(reflected_wave, shift1, order=3, mode="nearest", cval=0)
            transmitted_wave_new = shift(transmitted_wave, shift2, order=3, mode="nearest", cval=0)
            incident_end_new = -1*young_modulus*bar_area*(incident_wave + reflected_wave_new)
            transmission_end_new = -1*young_modulus*bar_area*(transmitted_wave_new)
            average_force_new = (incident_end_new+transmission_end_new)/2
            axs[0].clear()
            axs[0].plot(data['time'].values[:time_window_end]*1e6, incident_end_new, label='Incident End', color='b')
            axs[0].plot(data['time'].values[:time_window_end]*1e6, transmission_end_new, label='Transmission End', color='r')
            axs[0].plot(data['time'].values[:time_window_end]*1e6, average_force_new, label='Average', color='g')
            axs[0].set_xlabel('Time (us)')
            axs[0].set_ylabel('Force (N)')
            axs[0].set_title('SHPB Analysis')
            axs[0].legend()
            peak_index = reflected_wave_new.argmax()
            peak_time = data['time'].values[peak_index] * 1e6
            axs[1].clear()
            axs[1].plot(data['time'].values[:time_window_end]*1e6, incident_wave, label='Incident Wave', color='b')
            axs[1].plot(data['time'].values[:time_window_end]*1e6, reflected_wave_new, label='Reflected Wave',color='g')
            axs[1].axvline(peak_time, linestyle='--', color='grey')
            axs[1].plot(data['time'].values[:time_window_end]*1e6, transmitted_wave_new, label='Transmitted Wave',color='r')
            axs[1].set_xlabel('Time (us)')
            axs[1].set_ylabel('Amplitude')
            axs[1].set_title('Waves')
            axs[1].legend()
            fig.canvas.draw_idle()

        shift1_slider.on_changed(update)
        shift2_slider.on_changed(update)
        plt.show()

        #assign updated shifts

        if shift1 == optimal_shift1: 
            pass
        else:
            reflected_wave = shift(reflected_wave, shift1, order=3, mode="nearest", cval=0)

        if shift2 == optimal_shift2:
            pass
        else:
            transmitted_wave = shift(transmitted_wave, shift2, order=3, mode="nearest", cval=0)

        incident_end = -1*young_modulus*bar_area*(incident_wave + reflected_wave)
        transmission_end = -1*young_modulus*bar_area*(transmitted_wave)
        average_force = (incident_end+transmission_end)/2

        #####################################
        #analysis
        #####################################


        #print wave speed
        wave_speed = np.sqrt(young_modulus/density)

        #sample density
        sample_density = sample_mass/(sample_area*sample_thickness)

        #one wave analysis
        #strain
        sample_strain = 2*(wave_speed/sample_thickness)*cumtrapz(reflected_wave, data["time"].values[:time_window_end])

        #strain rate
        sample_strain_rate = 2*(wave_speed/sample_thickness)*reflected_wave

        strain_rate_argmax = sample_strain_rate.argmax()

        for i in range(len(sample_strain_rate)):
            if i > strain_rate_argmax:
                if sample_strain_rate[i] < 0:
                    strain_rate_cutoff = i-1
                    break
                else:
                    strain_rate_cutoff = len(sample_strain_rate)-1

        #stress strain curve

        sample_stress = -1*young_modulus*(bar_area/sample_area)*transmitted_wave
        #print(sample_stress)

        ultimate_tensile_strength = max(sample_stress)
        five_percent_uts = 0.05*ultimate_tensile_strength
        stress_argmax = sample_stress.argmax()
        strain_argmax = sample_strain.argmax()

        for i in range(len(sample_strain)-1):
            if sample_strain[i+1] < sample_strain[i]:
                if i > strain_argmax:
                    cut_off_index = i
                    break
            if sample_stress[i] < five_percent_uts:
                if stress_argmax < i:
                    cut_off_index = i
                    break

        #print(cut_off_index)
        
        #strain energy, specific energy

        strain_energy = cumtrapz(sample_stress[:cut_off_index], sample_strain[:cut_off_index])

        specific_energy = strain_energy/sample_density

        #####################################
        #REPORT GENERATION
        #####################################

        # Create a dictionary of dataframes, one for each plot
        data_dict = {
            'Force vs Time': pd.DataFrame({'time': data['time'].values[:time_window_end]*1e6, 'incident_end': incident_end, 'transmission_end': transmission_end, 'average_force': average_force}),
            'Strain vs Time': pd.DataFrame({'time': data['time'].values[:strain_rate_cutoff]*1e6, 'strain': sample_strain[:strain_rate_cutoff]}),
            'Strain Rate vs Time': pd.DataFrame({'time': data['time'].values[:strain_rate_cutoff]*1e6, 'strain_rate': sample_strain_rate[:strain_rate_cutoff]}),
            'Stress vs Strain': pd.DataFrame({'strain': sample_strain[:cut_off_index], 'stress': sample_stress[:cut_off_index]}),
            'Specific Energy vs Strain': pd.DataFrame({'strain': sample_strain[:cut_off_index-1], 'spec_energy': specific_energy}),
            'Strain Energy Density vs Strain': pd.DataFrame({'time': sample_strain[:cut_off_index-1], 'energy_density': strain_energy})
        }

        def export_to_excel(event):
            base_file_name = os.path.basename(file_path)
            base_file_name_without_extension = os.path.splitext(base_file_name)[0]
            excel_file_name = f"{base_file_name_without_extension}_excel.xlsx"
            writer = pd.ExcelWriter(excel_file_name, engine='xlsxwriter')
            # create a single sheet
            sheet_name = "Results"
            workbook = writer.book
            worksheet = workbook.add_worksheet(sheet_name)
            # write the data to the sheet, leaving a blank column to the right
            start_col = 0
            for section_name, df in data_dict.items():
                # Write the section name to the first row of the section
                worksheet.write(0, start_col, section_name)
                # Write the data to the sheet, leaving a blank column to the right
                df.to_excel(writer, sheet_name=sheet_name, startrow=1, startcol=start_col, index=False)
                # Increment the start column for the next section
                start_col += len(df.columns) + 1 # Add one for the blank column to the right
                # Insert a blank column to the right of the section
                worksheet.write_blank(0, start_col-1, "")

            # Set the column widths to be equal for all columns
            for i, col_width in enumerate([20] * start_col):
                worksheet.set_column(i, i, col_width)
                # Add a chart for the 'Stress vs Strain' data set
            writer.close()


        fig, axs = plt.subplots(2, 3, figsize=(10, 8))
        fig.subplots_adjust(hspace=0.4, wspace=0.4)
        axs[0, 0].plot(data['time'].values[:time_window_end]*1e6, incident_end, label='Incident End')
        axs[0, 0].plot(data['time'].values[:time_window_end]*1e6, transmission_end, label='Transmission End')
        axs[0, 0].plot(data['time'].values[:time_window_end]*1e6, average_force, label='Average')
        axs[0, 0].set(xlabel='Time (us)', ylabel='Force (N)', title='Stress Equilibrium')
        axs[0, 0].legend()

        axs[0, 1].plot(data['time'].values[:strain_rate_cutoff]*1e6, sample_strain[:strain_rate_cutoff])
        axs[0, 1].set(xlabel='Time (us)', ylabel='Strain', title='One Wave Strain vs. Time')

        axs[1, 0].plot(data['time'].values[:strain_rate_cutoff]*1e6, sample_strain_rate[:strain_rate_cutoff])
        axs[1, 0].set(xlabel='Time (us)', ylabel='Strain Rate', title='One Wave Strain Rate vs. Time')

        axs[1, 1].plot(sample_strain[:cut_off_index], sample_stress[:cut_off_index]/1e6)
        axs[1, 1].set(xlabel='Strain', ylabel='Stress (MPa)', title='Stress vs. Strain')

        axs[0, 2].plot(sample_strain[:cut_off_index-1], specific_energy/1e3)
        axs[0, 2].set(xlabel='Strain', ylabel='Specific Energy (kJ/kg)', title='Specific Energy vs. Strain')

        axs[1, 2].plot(sample_strain[:cut_off_index-1], strain_energy/1e3)
        axs[1, 2].set(xlabel='Strain', ylabel='Strain Energy Density (kJ/m^3)', title='Energy Density vs. Strain')

        ax_button = plt.axes([0.8, 0.02, 0.1, 0.035])
        button = Button(ax_button, 'Export to Excel')
        button.on_clicked(export_to_excel)

        plt.show()


