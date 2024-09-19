import os
import pandas as pd
import numpy as np
import plotly.graph_objs as go
from scipy.signal import find_peaks, savgol_filter, peak_widths
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Function to process the CSV file
def process_csv(file_path):
    # Read the original CSV file
    df = pd.read_csv(file_path, header=None)
    
    # Assign new column names for wave data
    wave_columns = ['wave0'] + [f'wave{i}' for i in range(1, df.shape[1])]

    # Extract temp1, temp2, and tim rows
    temp1 = df.iloc[0].values
    temp2 = df.iloc[1].values
    tim = df.iloc[2].values
    
    # Drop the first three rows from the original data (since they are temp1, temp2, tim)
    wave_data = df.iloc[3:].reset_index(drop=True)
    wave_data.columns = wave_columns
    
    # Create a new DataFrame with temp1, temp2, and tim as columns
    refactored_df = pd.DataFrame({
        'temp1': temp1,
        'temp2': temp2,
        'tim': tim
    })
    
    # Concatenate the temp1, temp2, tim columns with the wave data
    refactored_df = pd.concat([refactored_df, wave_data], axis=1)
    
    return refactored_df, wave_columns

# Function to save the refactored CSV
def save_refactored_csv(refactored_df, file_name, folder_path):
    # Create a folder for the original CSV file if it doesn't exist
    output_folder = os.path.join(folder_path, file_name)
    os.makedirs(output_folder, exist_ok=True)
    
    # Save the refactored CSV file inside the created folder
    refactored_file_path = os.path.join(output_folder, f"refactored_{file_name}.csv")
    refactored_df.to_csv(refactored_file_path, index=False)
    
    return refactored_file_path

# Function to find peaks and their widths
def find_peaks_and_widths(y_data, x_data):
    # Apply smoothing to avoid detecting peaks in non-smooth regions
    smoothed = savgol_filter(y_data, window_length=11, polyorder=3)
    
    # Find peaks
    peaks, properties = find_peaks(smoothed, prominence=0.0003, height=0, distance=1, width=2)
    
    # Calculate peak widths using the rel_height parameter (same approach as provided)
    widths_raw, width_heights, left_ips, right_ips = peak_widths(smoothed, peaks, rel_height=0.5)
    
    # Interpolate to get more precise x-values for start and end points
    start_x = np.interp(left_ips, np.arange(len(x_data)), x_data)
    end_x = np.interp(right_ips, np.arange(len(x_data)), x_data)
    
    # Calculate the width as the difference between start and end points
    widths = end_x - start_x
    
    return peaks, widths, start_x, end_x, width_heights

# Function to convert wavelength to point (pnt) in Igor style
def wavelength_to_pnt(wavelength, min_wavelength, max_wavelength, total_points):
    point_index = ( (wavelength - min_wavelength) / (max_wavelength - min_wavelength) ) * (total_points)
    return int(round(point_index))

# Function to save peaks and widths to CSV and include point-based columns
def save_peaks_and_widths(x_data, y_data, peaks, widths, start_x, end_x, wave, output_folder, min_wavelength, max_wavelength, total_points):
    # Convert start and end wavelengths to points
    start_pnt = [wavelength_to_pnt(start, min_wavelength, max_wavelength, total_points) for start in start_x]
    end_pnt = [wavelength_to_pnt(end, min_wavelength, max_wavelength, total_points) for end in end_x]
    
    # Calculate pnt_width as the difference between end_pnt and start_pnt
    pnt_width = np.array(end_pnt) - np.array(start_pnt)
    
    # Create a DataFrame to store the peak properties and widths
    df_peaks = pd.DataFrame({
        'Type': ['Peak'] * len(peaks),
        'Position': x_data[peaks],
        'Start': start_x,
        'End': end_x,
        'Width': widths,  # Use the correct width calculation (end - start)
        'pnt_x': start_pnt,
        'pnt_y': end_pnt,
        'pnt_width': pnt_width
    })
    
    # Save the widths and peak data to a CSV file
    output_file = os.path.join(output_folder, f"{wave}_peaks_widths.csv")
    df_peaks.to_csv(output_file, index=False)
    print(f"Peak widths for {wave} saved as {output_file}")

# Function to generate plots
def generate_plot(df, wave_columns, file_name, folder_path, min_wavelength, max_wavelength, total_points):
    wave0 = df['wave0']  # Use wave0 as the x-axis
    
    # Create a Plotly figure
    fig = go.Figure()

    # Create a list of dropdown options for each wave (excluding wave0)
    dropdown_buttons = []
    
    # Loop through all wave columns except wave0
    for i, wave in enumerate(wave_columns[1:]):
        y_data = df[wave].values
        
        # Find peaks and widths
        peaks, widths, start_x, end_x, width_heights = find_peaks_and_widths(y_data, wave0)
        
        # Add the wave plot (initially set to invisible)
        fig.add_trace(go.Scatter(
            x=wave0, y=df[wave], mode='lines', visible=False, name=f'wave0 vs {wave}'
        ))
        
        # Add peaks (initially set to invisible)
        fig.add_trace(go.Scatter(
            x=wave0[peaks], y=df[wave].values[peaks], mode='markers', visible=False, name=f'{wave} Peaks', marker=dict(color='red', symbol='x')
        ))
        
        # Create a button for each wave to toggle visibility
        visibility = [False] * len(fig.data)
        visibility[2 * i] = True  # Make the selected wave visible
        visibility[2 * i + 1] = True  # Make the corresponding peaks visible
        
        dropdown_buttons.append({
            'label': f'wave0 vs {wave}',
            'method': 'update',
            'args': [{'visible': visibility},  # Set only the selected wave visible
                     {'title': f'Wave0 vs {wave}'}],  # Update plot title
        })

        # Save peak widths to a CSV file, and include the point-based columns
        save_peaks_and_widths(wave0.values, df[wave].values, peaks, widths, start_x, end_x, wave, os.path.join(folder_path, file_name), min_wavelength, max_wavelength, total_points)

    # Set the first wave as visible by default
    fig.data[0].visible = True
    fig.data[1].visible = True
    
    # Add the dropdown menu
    fig.update_layout(
        updatemenus=[{
            'buttons': dropdown_buttons,
            'direction': 'down',  # Dropdown menu direction
            'showactive': True,
            'x': 0.17,  # Position of dropdown
            'y': 1.15
        }]
    )
    
    # Update layout for better visualization
    fig.update_layout(
        title=f'Wave0 vs Other Waves for {file_name}',
        xaxis_title='Wave0',
        yaxis_title='Wave Values',
        showlegend=True,
        hovermode='closest'
    )
    
    # Save as an HTML file inside the corresponding folder
    output_file = os.path.join(folder_path, file_name, f"{file_name}_plot.html")
    fig.write_html(output_file)
    print(f"Plot saved as {output_file}")

# Main function to process all files in the folder
def process_files_in_folder(base_folder, min_wavelength, max_wavelength, total_points):
    # Define the folder path for 'Physics Research'
    folder_path = os.path.join(os.path.expanduser('~'), 'Documents', base_folder)
    
    # Loop through all CSV files in the folder
    for filename in os.listdir(folder_path):
        if filename.endswith('.csv'):
            file_path = os.path.join(folder_path, filename)
            file_name = filename.replace('.csv', '')
            print(f"Processing {file_path}")
            
            # Step 1: Process and refactor the original CSV
            refactored_df, wave_columns = process_csv(file_path)
            
            # Step 2: Save the refactored CSV in a new folder named after the original file
            refactored_file_path = save_refactored_csv(refactored_df, file_name, folder_path)
            
            # Step 3: Read the refactored CSV and generate the graph (with peak detection)
            refactored_df = pd.read_csv(refactored_file_path)
            generate_plot(refactored_df, wave_columns, file_name, folder_path, min_wavelength, max_wavelength, total_points)

# Run the function to process and refactor all CSVs in the folder
process_files_in_folder('Physics Research', 1475, 1575, 1000)

# Define the Gaussian function for curve fitting
def gauss(x, a, x0, sigma):
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

# Main function equivalent to trackgauss4
def trackgauss4(wpre, stnum, increment, endnum, wlen, suffix, start, stop, width, folder_path, original_file_name):
    imax = int(np.floor((endnum - stnum) / increment) + 1)
    
    # Initialize arrays to store peak values, errors, temperatures, and times
    peakuse = np.zeros(imax)
    erruse = np.zeros(imax)
    tempuse = np.zeros(imax)
    timeuse = np.zeros(imax)
    
    # Load the refactored CSV from the correct subfolder
    refactored_csv_path = os.path.join(folder_path, original_file_name, f'refactored_{original_file_name}.csv')
    df = pd.read_csv(refactored_csv_path)
    
    # Extract temperature and timing from CSV
    temperature = df['temp1'].values
    timing = df['tim'].values
    
    # Iterate over each wave (from stnum to endnum with given increment)
    for i in range(imax):
        wave_col = f"{wpre}{i * increment + stnum}"
        
        # Extract the relevant wave data
        tempwavey = df[wave_col]
        wlen_data = df[wlen]
        
        # Determine max location in the specified range (start to stop)
        maxloc = tempwavey[start:stop].idxmax()
        
        # Calculate the bounds for the left and right minimum search
        left_min_start = max(0, maxloc - int(1.5 * width))
        left_min_end = max(0, maxloc - int(0.5 * width))
        right_min_start = min(len(tempwavey), maxloc + int(0.5 * width))
        right_min_end = min(len(tempwavey), maxloc + int(1.5 * width))
        
        # Check if the slices are valid and not empty
        if left_min_start < left_min_end and right_min_start < right_min_end:
            leftminloc = tempwavey[left_min_start:left_min_end].idxmin()
            yleft = tempwavey[leftminloc]

            rightminloc = tempwavey[right_min_start:right_min_end].idxmin()
            yright = tempwavey[rightminloc]
            
            # Calculate slope and intercept for normalization
            slope = (yright - yleft) / (rightminloc - leftminloc)
            intcpt = yright - slope * rightminloc
            
            # Normalize the wave data
            tempwaveynorm = tempwavey / (slope * wlen_data + intcpt)
            
            # Gaussian fit: Use a better guess for initial parameters
            try:
                # Extract the portion of the data around the peak for fitting
                xdata = wlen_data[maxloc-int(width/2):maxloc+int(width/2)]
                ydata = tempwaveynorm[maxloc-int(width/2):maxloc+int(width/2)]
                
                # Provide initial guess for [amplitude, mean, std deviation]
                initial_guess = [np.max(ydata), np.mean(xdata), np.std(xdata)]
                
                # Perform the Gaussian curve fit
                popt, pcov = curve_fit(gauss, xdata, ydata, p0=initial_guess)
                
                # Store peak, error, temperature, and timing
                peakuse[i] = popt[1]  # x0 from Gaussian fit
                erruse[i] = np.sqrt(np.diag(pcov))[1]  # error in x0 (std deviation)
                tempuse[i] = temperature[i * increment + stnum - 1]
                timeuse[i] = timing[i * increment + stnum - 1]
                
            except RuntimeError as e:
                print(f"Fit failed for wave {wave_col}: {e}")
        else:
            print(f"Warning: Invalid slice for wave {wave_col}")
    
    # Create the interactive Plotly graph with dual y-axes
    fig = go.Figure()

    # Add the peak wavelength (nm) on the left y-axis
    fig.add_trace(go.Scatter(
        x=timeuse, y=peakuse, mode='lines', name=f'Peak {suffix} (Wavelength)',
        line=dict(color='red')
    ))

    # Add the temperature (°C) on the right y-axis
    fig.add_trace(go.Scatter(
        x=timeuse, y=tempuse, mode='lines', name=f'Temperature {suffix} (°C)',
        line=dict(color='blue'), yaxis='y2'
    ))

    # Update layout for dual y-axes
    fig.update_layout(
        title=f'Peak Wavelength and Temperature for Peak {suffix}',
        xaxis_title='Time (seconds)',
        yaxis=dict(
            title='Peak Wavelength (nm)',
            titlefont=dict(color='red'),
            tickfont=dict(color='red')
        ),
        yaxis2=dict(
            title='Temperature (°C)',
            titlefont=dict(color='blue'),
            tickfont=dict(color='blue'),
            overlaying='y',
            side='right'
        ),
        legend=dict(x=0.01, y=0.99)
    )

    # Save the plot as an interactive HTML file
    output_html_path = os.path.join(folder_path, original_file_name, "gaussian fit", f"plot_{suffix}.html")
    os.makedirs(os.path.dirname(output_html_path), exist_ok=True)
    fig.write_html(output_html_path)
    print(f"Interactive plot saved as {output_html_path}")

# Process all peaks for each CSV file in the folder
def process_peaks_for_each_csv(folder_path, original_file_name, wpre="wave", wlen="wave0"):
    # Path to the folder containing the peak widths CSVs
    csv_folder = os.path.join(folder_path, original_file_name)
    
    # Iterate over all peak widths CSV files
    for peak_csv in os.listdir(csv_folder):
        if peak_csv.endswith("_peaks_widths.csv"):
            # Load the peak widths CSV
            peak_suffix = peak_csv.replace("_peaks_widths.csv", "")
            peak_df = pd.read_csv(os.path.join(csv_folder, peak_csv))
            
            # Extract pnt values for each peak
            for _, row in peak_df.iterrows():
                suffix = str(int(row['Position']))  # Use the peak position as the suffix
                start = int(row['pnt_x'])
                stop = int(row['pnt_y'])
                width = int(row['pnt_width'])
                
                # Call trackgauss4 for each peak
                trackgauss4(wpre, 1, 1, 546, wlen, suffix, start, stop, width, folder_path, original_file_name)

# Function to process all CSV files in the 'Physics Research' folder
def process_all_csv_files_in_folder(base_folder):
    folder_path = os.path.join(os.path.expanduser('~'), 'Documents', base_folder)
    
    # Iterate through each CSV file in the 'Physics Research' folder
    for filename in os.listdir(folder_path):
        if filename.endswith('.csv'):
            file_name = filename.replace('.csv', '')
            print(f"Processing file: {file_name}")
            
            # Create a new folder for each CSV file
            csv_folder = os.path.join(folder_path, file_name)
            os.makedirs(csv_folder, exist_ok=True)
            
            # Process all peaks for each CSV
            process_peaks_for_each_csv(folder_path, file_name)

# Run the function to process all CSV files in the 'Physics Research' folder
process_all_csv_files_in_folder('Physics Research')


