import os, glob, pyexcel
import matplotlib.pyplot as plt
from numpy import *
from scipy.optimize import curve_fit
from scipy.signal import medfilt

# Read  Me
# Specify the main directory below under workDir
# Then place data in the specific order inside the main directory : Sensor/Scenario/execfile
# For example: Strain_Gauge/SG_dyn/dv1.xlsx

# Main working directory
workDir = 'c:/Users/Richard/Desktop/Data_Reader'
# Change to main folder
os.chdir(workDir)
print('Now in {}'.format(workDir))

# Check all valid sensor folder
senList = [sen for sen in os.listdir(workDir) if os.path.isdir(sen)]
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

        # Declaring main variables arrays
        totalDisplacement = empty([3,50000])
        voltage = empty([3,50000])
        relativeDisplacement = empty([3,50000])

        # Iterating through 3 trials
        for trial in range(1,4):
            print('Reading trial {}'.format(trial))
            trialList = ['td'+str(trial)+'.xlsx', 'ta'+str(trial)+'.xlsx',
                        'dv'+str(trial)+'.xlsx',]
            # Total Displacement from td.xlsx
            # Time from ta.xlsx
            # Voltage from ta.xlsx
            # Relative displacement from dv.xlsx

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
                    totalDisplacement[trial-1] = data[1]*25.4
                    time = data[0]

                elif fileName == trialList[1]:
                    voltage[trial-1] = medfilt(data[1])


                elif  fileName == trialList[2]:
                    relativeDisplacement[trial-1] = data[0]*25.4

                else:
                    raise ValueError('Error reading the wrong trial file')

            # Filter all data
            def filtering():

                pass

            def rezeroing():
                global totalDisplacement, relativeDisplacement
                shift = average(totalDisplacement[trial-1])
                totalDisplacement[trial-1]-=shift
                relativeDisplacement[trial-1]-=shift

            # Curve fitting based on scenario
            def checkScenario(f, x, y, t):
                global report, trial, sensor

                if 'static' in folderName:
                    popt, pcov = curve_fit(f, x[trial-1], y[trial-1])
                    report[sensor]['static'].append(popt)

                elif 'qui' in folderName:
                    popt, pcov = curve_fit(f, x[trial-1], y[trial-1])
                    report[sensor]['quasi_static'].append(popt)

                elif 'tran' in folderName:
                    # print(len(x[trial-1][t>=1]))
                    popt, pcov = curve_fit(f, x[trial-1][t>=1], y[trial-1][t>=1])
                    # popt, pcov = curve_fit(f, x[trial-1], y[trial-1])
                    report[sensor]['transient'].append(popt)

                elif 'dyn' in folderName:
                    rezeroing()
                    popt, pcov = curve_fit(f, x[trial-1], y[trial-1])
                    report[sensor]['dynamic'].append(popt)

                else:
                    raise ValueError('Scenario Not Found ')

            # Data Anaysis
            if sensor == 'Strain_Gauge':
                # Governing equation
                def sg(x, a, b):
                    # Vout = Vin + g*IR*delta_L + I*R
                    return 5 + 2*a*x + b

                checkScenario(sg, relativeDisplacement, voltage, time)

            elif sensor == 'IR':
                # Governing equation
                def ir(x, a, b):
                    return a*(x**-2) + b

                checkScenario(ir, totalDisplacement, voltage, time)

            elif sensor == 'Voice_Coil':
                pass
            elif sensor == 'Accelerometer':
                pass
            else:
                raise ValueError('Sensor Not Regonized')

            # Plotting
            plt.figure(trial, figsize = (13, 7.5))
            plt.subplots(3,1)

            plt.subplot(311)
            plt.gca().set_title('Sensor: {} Case : {} Trial {}'.format(
            sensor, folderName, str(trial)))
            plt.plot(time, voltage[trial-1], label = '{}'.format(sensor))
            plt.xlabel('Time')
            plt.ylabel('Voltage')

            plt.subplot(312)
            plt.plot(time, totalDisplacement[trial-1], label = 'Encoder')
            plt.xlabel('Time')
            plt.ylabel('Distance')

            plt.subplot(313)
            plt.plot(totalDisplacement[trial-1], voltage[trial-1], label = 'vs')
            plt.xlabel('Distance')
            plt.ylabel('Voltage')

            plt.savefig('Trial_{}'.format(trial))
            plt.close()

    # Calculating Means
    # if sensor == 'Strain_Gauge':
    #     report[sensor]['static'] = append(report[sensor]['static'], average(report[sensor]['static']))
    #     report[sensor]['quasi_static'] = append(report[sensor]['quasi_static'], average(report[sensor]['quasi_static']))
    #     report[sensor]['transient'] = append(report[sensor]['transient'], average(report[sensor]['transient']))
    #     report[sensor]['dynamic'] = append(report[sensor]['dynamic'], average(report[sensor]['dynamic']))
    print(report[sensor])
