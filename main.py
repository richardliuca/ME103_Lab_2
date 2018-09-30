import os, glob, pyexcel
import matplotlib.pyplot as plt
from numpy import *
from scipy.optimize import curve_fit
from scipy.signal import medfilt, butter, filtfilt, freqz, get_window

# Read Me
# Specify the main directory below under workDir
# Then place data in the specific order inside the main directory : Sensor/Scenario/execfile
# For example: Strain_Gauge/SG_dyn/dv1.xlsx

# Governing equation
# Vout = Vin + g*IR*delta_L + I*R
sg = lambda x, a, b: 5 + 2*a*x + b
ir = lambda x, a, b: a*x + b
vc = lambda x, a, b: a*x + b
acc = lambda x, a, b: a*x + b

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
                'dynamic': [], 'Static_Model': [], 'Dynamic_Model': [],
                'xDataS': [], 'xDataD': [], 'xDataT': [],
                'yDataS': [], 'yDataD': [], 'yDataT': [],
                'residualS' : [], 'residualTD' : [], 'residualSD' : [],
                'StErrS' : [], 'StErrD': [],'SS_S' : [], 'SS_D': [],
                '95%CIS': [], '95%CID': []}
    # All Sensor summary list
    report[sensor] = sensDict

    # Check all valid scenrio folders
    dirList = [dir for dir in os.listdir(senDir) if os.path.isdir(dir)]
    for folderName in dirList:
        # Change to scenario folders
        os.chdir(senDir + '/' + folderName)
        print('Now in {}'.format(senDir+'/'+folderName))

        # residual = zeros((3, 50000))
        # xData = zeros((3, 50000))

        # Iterating through 3 trials
        for trial in range(1,4):
            # print('Reading trial {}'.format(trial))
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
                global time, voltage, totalDisplacement, sensor
                indexBool = time>=1
                voltage = voltage[indexBool]
                time = time[indexBool]
                if sensor == 'Accelerometer':
                    global acceleration
                    acceleration = acceleration[indexBool]
                elif sensor == 'Voice_Coil':
                    global velocity
                    velocity = velocity[indexBool]
                else:
                    totalDisplacement = totalDisplacement[indexBool]
                return indexBool

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
                global report, trial, sensor, residual, xData

                if 'static' in folderName:
                    y = low_pass_filter(y, cutoff = 10)
                    # y = average(x)*ones(len(x))
                    popt, _ = curve_fit(f, x, y)
                    report[sensor]['static'].append(popt)

                elif 'qui' in folderName:
                    y = low_pass_filter(y, cutoff = 20)
                    popt, _ = curve_fit(f, x, y)
                    report[sensor]['quasi_static'].append(popt)
                    report[sensor]['xDataS'].append(x)
                    report[sensor]['yDataS'].append(y)
                    report[sensor]['residualS'].append(y - f(x, popt[0], popt[1]))

                elif 'trans' in folderName:
                    index = truncate()
                    # if len(residual[0]) == 50000:
                    #     residual = zeros((3, len(x[index])))
                    #     xData = zeros((3, len(x[index])))
                    popt, _ = curve_fit(f, x[index], y[index])
                    report[sensor]['transient'].append(popt)
                    report[sensor]['xDataT'].append(x[index])
                    report[sensor]['yDataT'].append(y[index])
                    report[sensor]['residualTD'].append(y[index] - f(x[index], popt[0], popt[1]))

                elif 'dyn' in folderName:
                    rezeroing()
                    popt, _ = curve_fit(f, x, y)
                    report[sensor]['dynamic'].append(popt)
                    report[sensor]['xDataD'].append(x)
                    report[sensor]['yDataD'].append(y)
                    report[sensor]['residualSD'].append(y - f(x, popt[0], popt[1]))

                else:
                    raise ValueError('Scenario Not Found ')

                return popt

            # Plotting
            def plottingSubplot3(x1, y1, x2, y2, x3, y3, x4, y4,
                                xLabels, yLabels, labels):

                global sensor, folderName, trial
                plt.figure(figsize = (14, 7))

                plt.subplot(311)
                plt.gca().set_title('Sensor: {} Case : {} Trial {}'.format(
                sensor, folderName, str(trial)))
                plt.plot(x1, y1, label = labels[0])
                plt.xlabel(xLabels[0])
                plt.ylabel(yLabels[0])
                plt.grid(b = True, which = 'both')
                plt.legend()

                plt.subplot(312)
                plt.plot(x2, y2, label = labels[1])
                plt.xlabel(xLabels[1])
                plt.ylabel(yLabels[1])
                plt.grid(b = True, which = 'both')
                plt.legend()

                plt.subplot(313)
                plt.plot(x3, y3, label = labels[2])
                plt.plot(x4, y4, 'r', label = 'Fitted Curve')
                plt.xlabel(xLabels[2])
                plt.ylabel(yLabels[2])
                plt.grid(b = True, which = 'both')
                plt.legend()

                plt.savefig('Trial_{}.png'.format(trial))
                # plt.show()
                plt.close()

            # Data Anaysis
            if sensor == 'Strain_Gauge':
                coefficient = checkScenario(sg, totalDisplacement, voltage)

                fit = sg(totalDisplacement, coefficient[0], coefficient[1])

                # yBar = sum(voltage)/len(voltage)
                # ssreg = sum((fit - yBar)**2)
                # sstot = sum((voltage - yBar)**2)
                #
                # print('R square is : {}'.format(ssreg/sstot))

                # residual[trial-1] = voltage - fit
                # xData[trial-1] = totalDisplacement



                plottingSubplot3(time, voltage, time, totalDisplacement,
                            totalDisplacement, voltage,
                            totalDisplacement, fit,
                            ['Time (s)', 'Time (s)', 'Distance (m)'],
                            ['Voltage (v)', 'Distance (m)', 'Voltage (v)'],
                            [sensor, 'Encoder', 'data'])

            elif sensor == 'IR':
                coefficient = checkScenario(ir, totalDisplacement, voltage)

                fit = ir(totalDisplacement, coefficient[0], coefficient[1])

                # yBar = sum(voltage)/len(voltage)
                # ssreg = sum((fit - yBar)**2)
                # sstot = sum((voltage - yBar)**2)
                #
                # print('R square is : {}'.format(ssreg/sstot))
                #
                # residual[trial-1] = voltage - fit
                # xData[trial-1] = totalDisplacement

                plottingSubplot3(time, voltage, time, totalDisplacement,
                            totalDisplacement, voltage,
                            totalDisplacement, fit,
                            ['Time (s)', 'Time (s)', 'Distance (m)'],
                            ['Voltage (v)', 'Distance (m)', 'Voltage (v)'],
                            [sensor, 'Encoder', 'data'])

            elif sensor == 'Voice_Coil':
                velocity, _ = numerical_diff(totalDisplacement,
                                            dt = time[1]-time[0], scale = 40)
                coefficient = checkScenario(vc, velocity, voltage)

                fit = vc(velocity, coefficient[0], coefficient[1])

                # yBar = sum(voltage)/len(voltage)
                # ssreg = sum((fit - yBar)**2)
                # sstot = sum((voltage - yBar)**2)
                #
                # print('R square is : {}'.format(ssreg/sstot))
                #
                # residual[trial-1] = voltage -fit
                # xData[trial-1] = velocity

                plottingSubplot3(time, voltage, time, velocity,
                            velocity, voltage, velocity, fit,
                            ['Time (s)', 'Time (s)', 'Velocity (m/s)'],
                            ['Voltage (v)', 'Velocity (m/s)', 'Voltage (v)'],
                            [sensor, 'Encoder', 'data'])

            elif sensor == 'Accelerometer':
                velocity, _ = numerical_diff(totalDisplacement,
                                            dt = time[1]-time[0], scale = 40)
                acceleration, _ = numerical_diff(velocity,
                                    dt = time[1]-time[0], scale = 40)
                coefficient = checkScenario(acc, acceleration, voltage)

                fit = acc(acceleration, coefficient[0], coefficient[1])

                # yBar = sum(voltage)/len(voltage)
                # ssreg = sum((fit - yBar)**2)
                # sstot = sum((voltage - yBar)**2)
                #
                # print('R square is : {}'.format(ssreg/sstot))
                #
                    # residual[trial-1] = voltage - fit

                plottingSubplot3(time, voltage, time, acceleration,
                        acceleration, voltage, acceleration, fit,
                        ['Time (s)', 'Time (s)', 'Acceleration (m s^-2)'],
                        ['Voltage (v)', 'Acceleration (m s^-2)', 'Voltage (v)'],
                        [sensor, 'Encoder', 'data'])

            else:
                raise ValueError('Sensor Not Regonized')

    os.chdir(senDir)


    if sensor == 'Strain_Gauge' or sensor == 'IR' or sensor == 'Accelerometer':
        report[sensor]['Dynamic_Model'] = mean(report[sensor]['transient'] + report[sensor]['dynamic'], axis = 0)
        flatResD = append(reshape(report[sensor]['residualSD'], -1), reshape(report[sensor]['residualTD'], -1))
        flatXD = append(reshape(report[sensor]['xDataD'], -1) ,reshape(report[sensor]['xDataT'], -1))
        if sensor == 'Strain_Gauge' or sensor == 'IR':
            report[sensor]['Static_Model'] = mean(report[sensor]['quasi_static'], axis = 0)
            flatResS = reshape(report[sensor]['residualS'], -1)
            flatXS = reshape(report[sensor]['xDataS'], -1)
            report[sensor]['StErrS'] = sum(power(flatResS, 2)/(50000-2))
            report[sensor]['SS_S'] = sum(power(flatXS - mean(flatXS), 2))
            report[sensor]['95%CIS'] = [1.96*sqrt(report[sensor]['StErrS']/report[sensor]['SS_S']), sqrt(report[sensor]['StErrS']*(1/50000 + mean(flatXS)/report[sensor]['SS_S']))]

    elif sensor == 'Voice_Coil':
        report[sensor]['Dynamic_Model'] = mean(report[sensor]['transient'], axis = 0)
        flatResD = reshape(report[sensor]['residualTD'], -1)
        flatXD = reshape(report[sensor]['xDataT'], -1)

    report[sensor]['StErrD'] = sum(power(flatResD, 2)/(50000-2))
    report[sensor]['SS_D'] =  sum(power(flatXD - mean(flatXD), 2))
    report[sensor]['95%CID'] = [1.96*sqrt(report[sensor]['StErrD']/report[sensor]['SS_D']), sqrt(report[sensor]['StErrD']*(1/50000 + power(mean(flatXD), 2)/report[sensor]['SS_D']))]

    confBand = lambda x, StErr, SS, n, z: z*sqrt(StErr*(1 + 1/n + power(x - mean(x), 2)/SS))

    if sensor == 'Strain_Gauge':
        plt.figure(1, figsize = (14, 7))
        plt.title('Static Behavior : {}'.format(sensor))
        plt.plot(report[sensor]['xDataS'][0], report[sensor]['residualS'][0], 'k,', label = 'Trial 1 : Quasi-Transient')
        plt.plot(report[sensor]['xDataS'][1], report[sensor]['residualS'][1], 'k,', label = 'Trial 2 : Quasi-Transient')
        plt.plot(report[sensor]['xDataS'][2], report[sensor]['residualS'][2], 'k,', label = 'Trial 3 : Quasi-Transient')
        plt.xlabel('Independent (m)')
        plt.ylabel('Residual (v)')
        plt.grid(b = True, which = 'both')
        plt.legend()
        plt.savefig('Static_Behavior_Residual_{}.png'.format(sensor))
        # plt.show()
        plt.close()

        plt.figure(2, figsize = (14, 7))
        plt.title('Dynamic Behavior : {}'.format(sensor))
        plt.plot(report[sensor]['xDataT'][0], report[sensor]['residualTD'][0], 'k,', label = 'Trial 1 : Transient')
        plt.plot(report[sensor]['xDataT'][1], report[sensor]['residualTD'][1], 'k,', label = 'Trial 2 : Transient')
        plt.plot(report[sensor]['xDataT'][2], report[sensor]['residualTD'][2], 'k,', label = 'Trial 3 : Transient')

        plt.plot(report[sensor]['xDataD'][0], report[sensor]['residualSD'][0], 'k,', label = 'Trial 1 : Dynamic')
        plt.plot(report[sensor]['xDataD'][1], report[sensor]['residualSD'][1], 'k,', label = 'Trial 2 : Dynamic')
        plt.plot(report[sensor]['xDataD'][2], report[sensor]['residualSD'][2], 'k,', label = 'Trial 3 : Dynamic')
        plt.xlabel('Independent (m)')
        plt.ylabel('Residual (v)')
        plt.grid(b = True, which = 'both')
        plt.legend()
        plt.savefig('Dynamic_Behavior_Residual_{}.png'.format(sensor))
        # plt.show()
        plt.close()

        # print(confBand(flatXS, report[sensor]['StErrS'], report[sensor]['SS_S'], 50000, 1.96))
        print('Parameter Uncertanties Static: {} with 95% confidence'.format(report[sensor]['95%CIS']))
        plt.figure(3, figsize = (14, 7))
        plt.title('Static Model : {}'.format(sensor))
        plt.plot(flatXS, sg(flatXS, report[sensor]['Static_Model'][0], report[sensor]['Static_Model'][1]), 'k',label = 'Model')
        plt.plot(flatXS, sg(flatXS, report[sensor]['Static_Model'][0], report[sensor]['Static_Model'][1]) + confBand(flatXS, report[sensor]['StErrS'], report[sensor]['SS_S'], 50000, 1.96), 'r--',label = '95% Confidence Interval Upper Limit ')
        plt.plot(flatXS, sg(flatXS, report[sensor]['Static_Model'][0], report[sensor]['Static_Model'][1]) - confBand(flatXS, report[sensor]['StErrS'], report[sensor]['SS_S'], 50000, 1.96), 'r--',label = '95% Confidence Interval Lower Limit ')
        plt.plot(report[sensor]['xDataS'][0], report[sensor]['yDataS'][0], 'b,', label = 'Trial 1 : Quasi-Transient')
        plt.plot(report[sensor]['xDataS'][1], report[sensor]['yDataS'][1], 'b,', label = 'Trial 2 : Quasi-Transient')
        plt.plot(report[sensor]['xDataS'][2], report[sensor]['yDataS'][2], 'b,', label = 'Trial 3 : Quasi-Transient')
        plt.xlabel('Position (m)')
        plt.ylabel('Voltage (v)')
        plt.grid(b = True, which = 'both')
        plt.legend()
        plt.savefig('Static_Behavior_Model_{}.png'.format(sensor))
        # plt.show()
        plt.close()

        # print(confBand(flatXD, report[sensor]['StErrD'], report[sensor]['SS_D'], 50000, 1.96))
        print('Parameter Uncertanties Dynamic: {} with 95% confidence'.format(report[sensor]['95%CID']))
        plt.figure(4, figsize = (14, 7))
        plt.title('Dynamic Model : {}'.format(sensor))
        plt.plot(flatXD, sg(flatXD, report[sensor]['Dynamic_Model'][0], report[sensor]['Dynamic_Model'][1]), 'k',label = 'Model')
        plt.plot(flatXD, sg(flatXD, report[sensor]['Dynamic_Model'][0], report[sensor]['Dynamic_Model'][1]) + confBand(flatXD, report[sensor]['StErrD'], report[sensor]['SS_D'], 50000, 1.96), 'r--',label = '95% Confidence Interval Upper Limit ')
        plt.plot(flatXD, sg(flatXD, report[sensor]['Dynamic_Model'][0], report[sensor]['Dynamic_Model'][1]) - confBand(flatXD, report[sensor]['StErrD'], report[sensor]['SS_D'], 50000, 1.96), 'r--',label = '95% Confidence Interval Lower Limit ')
        plt.plot(report[sensor]['xDataT'][0], report[sensor]['yDataT'][0], 'b,', label = 'Trial 1 : Transient')
        plt.plot(report[sensor]['xDataT'][1], report[sensor]['yDataT'][1], 'b,', label = 'Trial 2 : Transient')
        plt.plot(report[sensor]['xDataT'][2], report[sensor]['yDataT'][2], 'b,', label = 'Trial 3 : Transient')

        plt.plot(report[sensor]['xDataD'][0], report[sensor]['yDataD'][0], 'b,', label = 'Trial 1 : Dynamic')
        plt.plot(report[sensor]['xDataD'][1], report[sensor]['yDataD'][1], 'b,', label = 'Trial 2 : Dynamic')
        plt.plot(report[sensor]['xDataD'][2], report[sensor]['yDataD'][2], 'b,', label = 'Trial 3 : Dynamic')
        plt.xlabel('Position (m)')
        plt.ylabel('Voltage (v)')
        plt.grid(b = True, which = 'both')
        plt.legend()
        plt.savefig('Dynamic_Behavior_Model_{}.png'.format(sensor))
        # plt.show()
        plt.close()

    elif sensor == 'IR':
        plt.figure(1, figsize = (14, 7))
        plt.title('Static Behavior : {}'.format(sensor))
        plt.plot(report[sensor]['xDataS'][0], report[sensor]['residualS'][0], 'k,', label = 'Trial 1 : Quasi-Transient')
        plt.plot(report[sensor]['xDataS'][1], report[sensor]['residualS'][1], 'k,', label = 'Trial 2 : Quasi-Transient')
        plt.plot(report[sensor]['xDataS'][2], report[sensor]['residualS'][2], 'k,', label = 'Trial 3 : Quasi-Transient')
        plt.xlabel('Independent (m)')
        plt.ylabel('Residual (v)')
        plt.grid(b = True, which = 'both')
        plt.legend()
        plt.savefig('Static_Behavior_Residual_{}.png'.format(sensor))
        # plt.show()
        plt.close()

        plt.figure(2, figsize = (14, 7))
        plt.title('Dynamic Behavior : {}'.format(sensor))
        plt.plot(report[sensor]['xDataT'][0], report[sensor]['residualTD'][0], 'k,', label = 'Trial 1 : Transient')
        plt.plot(report[sensor]['xDataT'][1], report[sensor]['residualTD'][1], 'k,', label = 'Trial 2 : Transient')
        plt.plot(report[sensor]['xDataT'][2], report[sensor]['residualTD'][2], 'k,', label = 'Trial 3 : Transient')

        plt.plot(report[sensor]['xDataD'][0], report[sensor]['residualSD'][0], 'k,', label = 'Trial 1 : Dynamic')
        plt.plot(report[sensor]['xDataD'][1], report[sensor]['residualSD'][1], 'k,', label = 'Trial 2 : Dynamic')
        plt.plot(report[sensor]['xDataD'][2], report[sensor]['residualSD'][2], 'k,', label = 'Trial 3 : Dynamic')
        plt.xlabel('Independent (m)')
        plt.ylabel('Residual (v)')
        plt.grid(b = True, which = 'both')
        plt.legend()
        plt.savefig('Dynamic_Behavior_Residual_{}.png'.format(sensor))
        # plt.show()
        plt.close()

        # print(confBand(flatXS, report[sensor]['StErrS'], report[sensor]['SS_S'], 50000, 1.96))
        print('Parameter Uncertanties Static: {} with 95% confidence'.format(report[sensor]['95%CIS']))
        plt.figure(3, figsize = (14, 7))
        plt.title('Static Model : {}'.format(sensor))
        plt.plot(flatXS, ir(flatXS, report[sensor]['Static_Model'][0], report[sensor]['Static_Model'][1]), 'k',label = 'Model')
        plt.plot(flatXS, ir(flatXS, report[sensor]['Static_Model'][0], report[sensor]['Static_Model'][1]) + confBand(flatXS, report[sensor]['StErrS'], report[sensor]['SS_S'], 50000, 1.96), 'r--',label = '95% Confidence Interval Upper Limit ')
        plt.plot(flatXS, ir(flatXS, report[sensor]['Static_Model'][0], report[sensor]['Static_Model'][1]) - confBand(flatXS, report[sensor]['StErrS'], report[sensor]['SS_S'], 50000, 1.96), 'r--',label = '95% Confidence Interval Lower Limit ')
        plt.plot(report[sensor]['xDataS'][0], report[sensor]['yDataS'][0], 'b,', label = 'Trial 1 : Quasi-Transient')
        plt.plot(report[sensor]['xDataS'][1], report[sensor]['yDataS'][1], 'b,', label = 'Trial 2 : Quasi-Transient')
        plt.plot(report[sensor]['xDataS'][2], report[sensor]['yDataS'][2], 'b,', label = 'Trial 3 : Quasi-Transient')
        plt.xlabel('Position (m)')
        plt.ylabel('Voltage (v)')
        plt.grid(b = True, which = 'both')
        plt.legend()
        plt.savefig('Static_Behavior_Model_{}.png'.format(sensor))
        # plt.show()
        plt.close()

        # print(confBand(flatXD, report[sensor]['StErrD'], report[sensor]['SS_D'], 50000, 1.96))
        print('Parameter Uncertanties Dynamic: {} with 95% confidence'.format(report[sensor]['95%CID']))
        plt.figure(4, figsize = (14, 7))
        plt.title('Dynamic Model : {}'.format(sensor))
        plt.plot(flatXD, ir(flatXD, report[sensor]['Dynamic_Model'][0], report[sensor]['Dynamic_Model'][1]), 'k',label = 'Model')
        plt.plot(flatXD, ir(flatXD, report[sensor]['Dynamic_Model'][0], report[sensor]['Dynamic_Model'][1]) + confBand(flatXD, report[sensor]['StErrD'], report[sensor]['SS_D'], 50000, 1.96), 'r--',label = '95% Confidence Interval Upper Limit ')
        plt.plot(flatXD, ir(flatXD, report[sensor]['Dynamic_Model'][0], report[sensor]['Dynamic_Model'][1]) - confBand(flatXD, report[sensor]['StErrD'], report[sensor]['SS_D'], 50000, 1.96), 'r--',label = '95% Confidence Interval Lower Limit ')
        plt.plot(report[sensor]['xDataT'][0], report[sensor]['yDataT'][0], 'b,', label = 'Trial 1 : Transient')
        plt.plot(report[sensor]['xDataT'][1], report[sensor]['yDataT'][1], 'b,', label = 'Trial 2 : Transient')
        plt.plot(report[sensor]['xDataT'][2], report[sensor]['yDataT'][2], 'b,', label = 'Trial 3 : Transient')

        plt.plot(report[sensor]['xDataD'][0], report[sensor]['yDataD'][0], 'b,', label = 'Trial 1 : Dynamic')
        plt.plot(report[sensor]['xDataD'][1], report[sensor]['yDataD'][1], 'b,', label = 'Trial 2 : Dynamic')
        plt.plot(report[sensor]['xDataD'][2], report[sensor]['yDataD'][2], 'b,', label = 'Trial 3 : Dynamic')
        plt.xlabel('Position (m)')
        plt.ylabel('Voltage (v)')
        plt.grid(b = True, which = 'both')
        plt.legend()
        plt.savefig('Dynamic_Behavior_Model_{}.png'.format(sensor))
        # plt.show()
        plt.close()

    elif sensor == 'Accelerometer':
        plt.figure(2, figsize = (14, 7))
        plt.title('Dynamic Behavior : {}'.format(sensor))
        plt.plot(report[sensor]['xDataT'][0], report[sensor]['residualTD'][0], 'k,', label = 'Trial 1 : Transient')
        plt.plot(report[sensor]['xDataT'][1], report[sensor]['residualTD'][1], 'k,', label = 'Trial 2 : Transient')
        plt.plot(report[sensor]['xDataT'][2], report[sensor]['residualTD'][2], 'k,', label = 'Trial 3 : Transient')

        plt.plot(report[sensor]['xDataD'][0], report[sensor]['residualSD'][0], 'k,', label = 'Trial 1 : Dynamic')
        plt.plot(report[sensor]['xDataD'][1], report[sensor]['residualSD'][1], 'k,', label = 'Trial 2 : Dynamic')
        plt.plot(report[sensor]['xDataD'][2], report[sensor]['residualSD'][2], 'k,', label = 'Trial 3 : Dynamic')
        plt.xlabel('Independent (m)')
        plt.ylabel('Residual (v)')
        plt.grid(b = True, which = 'both')
        plt.legend()
        plt.savefig('Dynamic_Behavior_Residual_{}.png'.format(sensor))
        # plt.show()
        plt.close()

        # print(confBand(flatXD, report[sensor]['StErrD'], report[sensor]['SS_D'], 50000, 1.96))
        print('Parameter Uncertanties Dynamic: {} with 95% confidence'.format(report[sensor]['95%CID']))
        plt.figure(3, figsize = (14, 7))
        plt.title('Dynamic Model : {}'.format(sensor))
        plt.plot(flatXD, acc(flatXD, report[sensor]['Dynamic_Model'][0], report[sensor]['Dynamic_Model'][1]), 'k',label = 'Model')
        plt.plot(flatXD, acc(flatXD, report[sensor]['Dynamic_Model'][0], report[sensor]['Dynamic_Model'][1]) + confBand(flatXD, report[sensor]['StErrD'], report[sensor]['SS_D'], 50000, 1.96), 'r--',label = '95% Confidence Interval Upper Limit ')
        plt.plot(flatXD, acc(flatXD, report[sensor]['Dynamic_Model'][0], report[sensor]['Dynamic_Model'][1]) - confBand(flatXD, report[sensor]['StErrD'], report[sensor]['SS_D'], 50000, 1.96), 'r--',label = '95% Confidence Interval Lower Limit ')
        plt.plot(report[sensor]['xDataT'][0], report[sensor]['yDataT'][0], 'b,', label = 'Trial 1 : Transient')
        plt.plot(report[sensor]['xDataT'][1], report[sensor]['yDataT'][1], 'b,', label = 'Trial 2 : Transient')
        plt.plot(report[sensor]['xDataT'][2], report[sensor]['yDataT'][2], 'b,', label = 'Trial 3 : Transient')

        plt.plot(report[sensor]['xDataD'][0], report[sensor]['yDataD'][0], 'b,', label = 'Trial 1 : Dynamic')
        plt.plot(report[sensor]['xDataD'][1], report[sensor]['yDataD'][1], 'b,', label = 'Trial 2 : Dynamic')
        plt.plot(report[sensor]['xDataD'][2], report[sensor]['yDataD'][2], 'b,', label = 'Trial 3 : Dynamic')
        plt.xlabel('Acceleration (m s^-2)')
        plt.ylabel('Voltage (v)')
        plt.grid(b = True, which = 'both')
        plt.legend()
        plt.savefig('Dynamic_Behavior_Model_{}.png'.format(sensor))
        # plt.show()
        plt.close()

    elif sensor == 'Voice_Coil':
        plt.figure(2, figsize = (14, 7))
        plt.title('Dynamic Behavior : {}'.format(sensor))
        plt.plot(report[sensor]['xDataT'][0], report[sensor]['residualTD'][0], 'k,', label = 'Trial 1 : Transient')
        plt.plot(report[sensor]['xDataT'][1], report[sensor]['residualTD'][1], 'k,', label = 'Trial 2 : Transient')
        plt.plot(report[sensor]['xDataT'][2], report[sensor]['residualTD'][2], 'k,', label = 'Trial 3 : Transient')
        plt.xlabel('Independent (m)')
        plt.ylabel('Residual (v)')
        plt.grid(b = True, which = 'both')
        plt.legend()
        plt.savefig('Dynamic_Behavior_Residual_{}.png'.format(sensor))
        # plt.show()
        plt.close()

        # print(confBand(flatXD, report[sensor]['StErrD'], report[sensor]['SS_D'], 50000, 1.96))
        print('Parameter Uncertanties Dynamic: {} with 95% confidence'.format(report[sensor]['95%CID']))
        plt.figure(3, figsize = (14, 7))
        plt.title('Dynamic Model : {}'.format(sensor))
        plt.plot(flatXD, vc(flatXD, report[sensor]['Dynamic_Model'][0], report[sensor]['Dynamic_Model'][1]), 'k',label = 'Model')
        plt.plot(flatXD, vc(flatXD, report[sensor]['Dynamic_Model'][0], report[sensor]['Dynamic_Model'][1]) + confBand(flatXD, report[sensor]['StErrD'], report[sensor]['SS_D'], 50000, 1.96), 'r--',label = '95% Confidence Interval Upper Limit ')
        plt.plot(flatXD, vc(flatXD, report[sensor]['Dynamic_Model'][0], report[sensor]['Dynamic_Model'][1]) - confBand(flatXD, report[sensor]['StErrD'], report[sensor]['SS_D'], 50000, 1.96), 'r--',label = '95% Confidence Interval Lower Limit ')
        plt.plot(report[sensor]['xDataT'][0], report[sensor]['yDataT'][0], 'b,', label = 'Trial 1 : Transient')
        plt.plot(report[sensor]['xDataT'][1], report[sensor]['yDataT'][1], 'b,', label = 'Trial 2 : Transient')
        plt.plot(report[sensor]['xDataT'][2], report[sensor]['yDataT'][2], 'b,', label = 'Trial 3 : Transient')
        plt.xlabel('Velocity (m/s)')
        plt.ylabel('Voltage (v)')
        plt.grid(b = True, which = 'both')
        plt.legend()
        plt.savefig('Dynamic_Behavior_Model_{}.png'.format(sensor))
        # plt.show()
        plt.close()

    else:
        pass

print('Complete')
