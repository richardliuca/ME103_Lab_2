import os, glob, pyexcel
import matplotlib.pyplot as plt
from numpy import *
from scipy.optimize import curve_fit
from scipy.signal import medfilt, butter, filtfilt, freqz, get_window

# Read  Me
# Specify the main directory below under workDir
# Then place data in the specific order inside the main directory : Sensor/Scenario/execfile
# For example: Strain_Gauge/SG_dyn/dv1.xlsx

# Main working directory
workDir = 'c:/Users/Richard/Documents/GitHub/ME103_Lab_2'
# Change to main folder
os.chdir(workDir)
print('Now in {}'.format(workDir))

# Check all valid sensor folder
senList = [sen for sen in os.listdir(workDir) if os.path.isdir(sen)]
senList.remove('.git')
senList.remove('__pycache__')
report = dict.fromkeys(senList)
for sensor in senList:
    senDir = workDir + '/' + sensor
    os.chdir(senDir)
    print('Now in {}'.format(senDir))
    # Sensor Summary Template
    sensDict = {'static': [], 'quasi_static': [], 'transient': [],
                'dynamic': [], 'mean': []}
    # All Sensor summary list
    report[sensor] = sensDict

    # Check all valid scenrio folders
    dirList = [dir for dir in os.listdir(senDir) if os.path.isdir(dir)]
    for folderName in dirList:
        # Change to scenario folders
        os.chdir(senDir + '/' + folderName)
        print('Now in {}'.format(senDir+'/'+folderName))

        # Iterating through 3 trials
        for trial in range(1,4):
            print('Reading trial {}'.format(trial))
            trialList = ['td'+str(trial)+'.xlsx', 'ta'+str(trial)+'.xlsx']
            # Total Displacement from td.xlsx
            # Time from ta.xlsx
            # Voltage from ta.xlsx

            # Low pass filter
            def low_pass_filter(x, cutoff = 124, sampling = 10000, order = 6):
                nyquist = sampling/2
                normalize_freq = (cutoff/nyquist)
                b, a = butter(order, normalize_freq,
                            btype = 'lowpass', analog = False)
                # Filter Frequency Response
                # w, h = freqz(b, a)
                # plt.semilogx(w/(2*pi)*sampling, 20*log10(abs(h)))
                # plt.axvline(12.4, color='k')
                # plt.grid(which='both', axis='both')
                # plt.title('Butterworth filter frequency response')
                # plt.xlabel('Frequency [Hz]')
                # plt.ylabel('Amplitude [dB]')
                # plt.show()


                x[:] = filtfilt(b, a, x)
                return x
                # return filtfilt(b, a, x)

            # Iterating throught data files to collect data
            for fileName in trialList:
                # Getting sheets from file
                sheet = pyexcel.get_sheet(file_name = fileName)
                # Converting data in sheet to array (50001x2)
                rawData = sheet.to_array()
                # Removing axis entries
                del rawData[0]
                # Transform to 2x50000 matrix
                data = array(rawData, dtype = float).transpose()
                # Selectively collecting data from each file
                if fileName == trialList[0]:
                    time = data[0]
                    totalDisplacement = data[1]*0.0254

                elif fileName == trialList[1]:
                    voltage = low_pass_filter(medfilt(medfilt(data[1])))

                else:
                    raise ValueError('Error reading the wrong trial file')

            def rezeroing():
                global totalDisplacement
                shift = average(totalDisplacement)
                totalDisplacement-=shift

            def truncate():
                global time, voltage, totalDisplacement
                indexBool = time>=1
                voltage = voltage[indexBool]
                totalDisplacement = totalDisplacement[indexBool]
                time = time[indexBool]

            def numerical_diff(x, dt, scale):
                global low_pass_filter
                xdot = zeros(len(x))
                xdotdot = zeros(len(x))
                for i in range(len(x)):
                    if i <= 1*scale:
                        # Forward differentiation
                        xdot[i] = (-3*x[i] + 4*x[i+1*scale]
                                    - x[i+2*scale])/(2*scale*dt)
                        xdotdot[i] = (x[i+2*scale] - 2*x[i+1*scale]
                                    + x[i])/((scale*dt)**2)
                    elif i >= (len(x) - 2*scale):
                        # Backward differentiation
                        xdot[i] = (3*x[i] - 4*x[i-1*scale]
                                    + x[i-2*scale])/(2*scale*dt)
                        xdotdot[i] = (-x[i-2*scale] + 2*x[i-1*scale]
                                    - x[i])/((scale*dt)**2)
                    else:
                        # Fourth order center differentiation
                        xdot[i] = (-x[i+2*scale]
                                    + 8*x[i+1*scale] - 8*x[i-1*scale]
                                    + x[i-2*scale])/(12*scale*dt)
                        xdotdot[i] = (-x[i+2*scale] + 16*x[i+1*scale] -30*x[i]
                                    + 16*x[i-1*scale]
                                    - x[i-2*scale])/(12*(scale*dt)**2)

                return low_pass_filter(medfilt(xdot)), xdotdot

            # Curve fitting based on scenario
            def checkScenario(f, x, y):
                global report, trial, sensor

                if 'static' in folderName:
                    y = low_pass_filter(y, cutoff = 10)
                    # y = average(x)*ones(len(x))
                    popt, _ = curve_fit(f, x, y)
                    report[sensor]['static'].append(popt)

                elif 'qui' in folderName:
                    y = low_pass_filter(y, cutoff = 20)
                    popt, _ = curve_fit(f, x, y)
                    report[sensor]['quasi_static'].append(popt)

                elif 'trans' in folderName:
                    truncate()
                    popt, _ = curve_fit(f, x, y)
                    report[sensor]['transient'].append(popt)

                elif 'dyn' in folderName:
                    rezeroing()
                    popt, _ = curve_fit(f, x, y)
                    report[sensor]['dynamic'].append(popt)

                else:
                    raise ValueError('Scenario Not Found ')

            # Data Anaysis
            if sensor == 'Strain_Gauge':
                # Governing equation
                # Vout = Vin + g*IR*delta_L + I*R

                sg = lambda x, a, b: 5 + 2*a*x + b
                checkScenario(sg, totalDisplacement, voltage)

                # sg = lambda x, a, b: a*x + b
                # checkScenario(sg, voltage, totalDisplacement)

            elif sensor == 'IR':
                # Governing equation
                ir = lambda x, a, b: a*exp(b*x)
                checkScenario(ir, totalDisplacement, voltage)

                # ir = lambda x, a, b: a*sqrt(x) + b
                # checkScenario(ir, voltage, totalDisplacement)

            elif sensor == 'Voice_Coil':
                pass
                # velocity, _ = numerical_diff(totalDisplacement,
                #                             dt = time[1]-time[0], scale = 40)
                # plt.figure(5)
                # plt.title('Encoder Velocity')
                # plt.plot(time, velocity, 'b')
                # plt.xlabel('Time (s)')
                # plt.ylabel('Velocity (mm/s)')
                # plt.show()

            elif sensor == 'Accelerometer':
                pass
                # velocity, _ = numerical_diff(totalDisplacement,
                #                             dt = time[1]-time[0], scale = 40)
                # acceleration, _ = numerical_diff(velocity,
                #                     dt = time[1]-time[0], scale = 40)
                # plt.figure(4)
                # plt.title('Encoder Acceleration')
                # plt.plot(time, acceleration, 'b')
                # plt.xlabel('Time (s)')
                # plt.ylabel('Acceleration (mm s^-2)')
                # plt.show()

            else:
                raise ValueError('Sensor Not Regonized')

            # Plotting
            plt.figure(figsize = (13, 7.5))
            # plt.subplots(3,1)

            plt.subplot(311)
            plt.gca().set_title('Sensor: {} Case : {} Trial {}'.format(
            sensor, folderName, str(trial)))
            plt.plot(time, voltage, label = '{}'.format(sensor))
            plt.xlabel('Time')
            plt.ylabel('Voltage')
            plt.legend()

            plt.subplot(312)
            plt.plot(time, totalDisplacement, label = 'Encoder')
            plt.xlabel('Time')
            plt.ylabel('Distance')
            plt.legend()

            plt.subplot(313)
            plt.plot(totalDisplacement, voltage, label = 'curve fit')
            # plt.plot(voltage, totalDisplacement, label = 'curve fit')
            # plt.xlabel('Voltage')
            plt.xlabel('Distance')
            # plt.ylabel('Distance')
            plt.ylabel('Voltage')
            plt.legend()

            plt.savefig('Trial_{}'.format(trial))
            # plt.show()
            plt.close()

    # Calculating Means
    # if sensor == 'Strain_Gauge':
    #     report[sensor]['static'] = append(report[sensor]['static'], average(report[sensor]['static']))
    #     report[sensor]['quasi_static'] = append(report[sensor]['quasi_static'], average(report[sensor]['quasi_static']))
    #     report[sensor]['transient'] = append(report[sensor]['transient'], average(report[sensor]['transient']))
    #     report[sensor]['dynamic'] = append(report[sensor]['dynamic'], average(report[sensor]['dynamic']))
    print(report[sensor])
