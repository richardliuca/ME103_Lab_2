import pyexcel
from scipy.signal import medfilt
from low_pass_filter import *
from numpy import *
import matplotlib.pyplot as plt


print('Expecting excel with time in 1st col and voltage in 2nd col')
filename = input('Please enter the excel file name: ') + '.xlsx'
sheet = pyexcel.get_sheet(file_name = filename)
rawData = sheet.to_array()

# Removing axis entries
del rawData[0]
# Transform to 2x50000 matrix
data = array(rawData, dtype = float).transpose()
time = data[0]
voltage = data[1]

voltageClean = low_pass_filter(medfilt(medfilt(voltage)), cutoff = 10000*2.5, sampling = 100000)

plt.plot(time, voltage, 'b', label = 'raw data with noise')
plt.plot(time, voltageClean, 'r', label = 'low pass filtered data')
plt.legend()
plt.show()

export = array([time, voltageClean])

pyexcel.isave_as(array = export.transpose(), dest_file_name = 'filtered_data_.xlsx')
print('Filtered data is exported as filtered_data_.xlsx')
print('First column is time and second column is voltage')
