#!/bin/env python3

import numpy as np

#========#
#========#

def read_data(filename):

    mydata1 = {}
    header_dict = {}
    flag_header = False
    idcounter = 0
    linecount = 0
    flag_id = False
    flag_output = False
    output = []

    # todo flag_id remove ?

    # Determine EOF
    with open(filename) as f:
        f.readlines()
        numEOF = f.tell()
    f.close()

    with open(filename) as f:
        for line in f:
            if line=='\n' or line.startswith('='):

                if line.startswith('===t') and flag_header==False:
                    header_list = line.split()
                    header_list[0] = 'time'
                    for i,val in enumerate(header_list):
                        mydata1.update({val: []})
                        header_dict.update({val: []})
                    flag_header=True

                # last dataset
                elif line.startswith('===EOF'):
                    #--------------------------------------------------------------#
                    flag_id = False
                    flag_output = False

                    # new dataset in output
                    output.append({})

                    for i, val in enumerate(header_list):
                        output[idcounter].update({val:[]})
                        for i in range(len(mydata1[val])):
                            output[idcounter][val].append(mydata1[val][i])

                        mydata1[val].clear()
                    idcounter += 1
                    #--------------------------------------------------------------#
                else:
                    continue

            elif line.startswith('0.0000'):

                if flag_output == False:
                    #--------------------------------------------------------------#
                    tmp = [float(n) for n in line.split()]
                    for i,val in enumerate(header_list):
                        mydata1[val].append(tmp[i])
                    #--------------------------------------------------------------#
                    flag_output = True

                elif flag_output == True:
                    #--------------------------------------------------------------#
                    # new dataset in output
                    output.append({})

                    for i, val in enumerate(header_list):
                        output[idcounter].update({val: []})
                        for i in range(len(mydata1[val])):
                            output[idcounter][val].append(mydata1[val][i])

                        mydata1[val].clear()

                    #--------------------------------------------------------------#
                    tmp = [float(n) for n in line.split()]
                    for i,val in enumerate(header_list):
                        mydata1[val].append(tmp[i])
                    #--------------------------------------------------------------#

                    idcounter += 1
                    #--------------------------------------------------------------#

                # Now reading starts

            # plain input to mydata
            elif flag_output == True:

                tmp = [float(n) for n in line.split()]

                for i,val in enumerate(header_list):
                    mydata1[val].append(tmp[i])


    return output

#========#
#========#

def normalize_dict(inquiry: dict):

    if not inquiry:
        print('>>> Warning! Dictionary was empty! <<<')
        return {}

    inquiry_new = {}

    try:
        for i,value in enumerate(dict(inquiry)):
            inquiry_new[value] = inquiry[value]
    except TypeError:
        error_string = 'Not a dictionary type class!'
        raise TypeError(error_string)


    for i,(valA,valB) in enumerate(inquiry_new.items()):
        if type(valB)!=float and type(valB)!=int:
            raise ValueError(valB,'is not a number')
        if float(valB) < 0:
            print ('Input is negative. They are ignored!')
            continue

    sum = 0
    for i,(valA,valB) in enumerate(inquiry_new.items()):
        if valB < 0:
            valB = 0
        sum += valB

    for i,(valA,valB) in enumerate(inquiry_new.items()):
        inquiry_new[valA] = valB/sum

    return inquiry_new

#========#
#========#

def func_inter(timex,data):
    time_array = data[0]
    func_array = data[1]
    funcX = np.interp(timex,time_array,func_array)
    return float(funcX)

#========#
#========#

def convX_mass_xpos(xval, dataset, wProxy, wCVol0, wC0):
    # dataset of [time mp x y z]
    # wProxy (ar)
    # wC0 (ar)
    # wCVol0 (only Vol)

    mCx0 = 0    # mass of particle at
    mCx  = 0

    wFC0  = wProxy['wFC']
    wVol0 = wProxy['wVol']
    wh2o0 = wProxy['wH2O']
    wAsh0 = wProxy['wAsh']
    mC0  = dataset['mp'][0] * wC0

    mStage0 = dataset['mp'][0]
    mStage1 = (1-wh2o0) * dataset['mp'][0]
    mStage2 = (1-wh2o0-wVol0) * dataset['mp'][0]
    mStage3 = (1-wh2o0-wVol0-wFC0) * dataset['mp'][0]


    for i in range(len(dataset['x'])):

        if xval == 0.0:
            mCx += mC0
            mCx0 += mC0

        if i == 0:
            continue

        x0  = dataset['x'][i-1] - xval
        x1  = dataset['x'][i] - xval

        if ( x1 > 0 and x0 < 0 ) or ( x1 < 0 and x0 > 0 ):

            if  dataset['mp'][i] > mStage1:
                mCx += mC0
                mCx0 += mC0
            elif dataset['mp'][i] <= mStage1 and dataset['mp'][i] > mStage2:
                # mStage1 is amount of water
                # only part is carbon of vol gas! (here fixed part)
                mCx += mC0 - (mStage1 - func_inter(xval, [dataset['x'], dataset['mp']])) * wCVol0
                mCx0 += mC0
            elif dataset['mp'][i] <= mStage2 and dataset['mp'][i] > mStage3:
                # 100 percent of the massloss comes from carbon
                mCx += mC0 - (mStage1-mStage2)*wCVol0 - (mStage2 - func_inter(xval, [dataset['x'], dataset['mp']]))
                mCx0 += mC0
            elif dataset['mp'][i] < mStage3:
                mCx += 0
                mCx0 += mC0

    return mCx, mCx0


#####################################################
## START
#####################################################

#-------------------------------
data_names = [
    "./dpm-data-inj1-mod.his",
    "./dpm-data-inj2-mod.his",
    "./dpm-data-inj3-mod.his",
    "./dpm-data-inj4-mod.his",
    "./dpm-data-inj5-mod.his",
]


#----------------------------
#-- OUTPUT ------------------
#----------------------------
zeroData1 = []
for i,name in enumerate(data_names):
    #--------------------
    # zeroData1 .. total sets
    # zeroData1[] .. set xy
    # zeroData1[][] .. particle zz from set xy; inj1 with 15 particles
    # zeroData1[][]['time'] .. data stored in dict syntax
    # zeroData1[][]['time'][] .. actual datapoint
    # dict: time mp x y z
    # one .. total set
    # two .. individual sets
    zeroData1.append(read_data(data_names[i]))

    #--------------------

#-- Get the mass of mp at position x
reactor_len = 2.20

#wp={'wFC':0.5428,'wVol':0.4527,'wAsh':0.0,'wH2O':0.0}
#wp={'wFC':0.48993128,'wVol':0.41266872,'wAsh':0.0974,'wH2O':0.0}
wp={'wFC': 0.4678843724, 'wVol': 0.3940986276, 'wAsh': 0.093017, 'wH2O': 0.045}

wUlit    = {'C': 0.724983256305}       # ar
# please use vol without tar ref because otherwise you need to calc both together!
wUltiVol = {'C': 0.6434820647419073}   # ar

# as received fraction
#wUlit['C'] =  wUlit['C'] * (1-wp['wAsh']) * (1-wp['wH2O'])

#wUltiVol['C'] = wUltiVol['C'] * (1-wp['wAsh']) * (1-wp['wH2O'])
#wp['wFC']  = wp['wFC']  * (1-wp['wH2O'])
#wp['wVol'] = wp['wVol'] * (1-wp['wH2O'])
#wp['wAsh'] = wp['wAsh'] * (1-wp['wH2O'])

# just the check
#wp_check = normalize_dict(wp)

#for i,val in enumerate(wp):
#    diff = wp[val] - wp_check[val]
#    if diff != 0:
#        print('Diff in fraction dict of the particle:',wp[val],wp_check[val])




xrange = np.linspace(0.0, reactor_len, num=200)
convX_dict = {'x':[],'mC':[],'mC0':[],'XC':[]}

#convX_dict = {'x':[],'mp':[],'mp0':[],'mChar0':[],'mC0':[],'Xp':[]}
mass_tmp   = 0
massChar0_tmp = 0
massC_tmp = 0
massC0_tmp = 0
massp0_tmp = 0
counter_tmp = 0


#dud = convX_mass_xpos(1.10,zeroData1[4][0],wp,wUltiVol['C'],wUlit['C'])
#test = (dud[1]-dud[0])/dud[1]
print('xval   mC   mC0 ')

mweight = np.array([0.1,4/15,4/15,4/15,0.1])
#mweight = np.array([1,1,1,1,1])
#mweight = mweight #* 5.555e-4 * wUlit['C']

for i,xval in enumerate(xrange):
    for j in range(len(zeroData1)):  # injA

        #mweightx = mweight[j] / (float(len(zeroData1[j])) * zeroData1[j][0]['mp'][0])
        mweightx = mweight[j] / zeroData1[j][0]['mp'][0]

        for k in range(len(zeroData1[j])):  # particle A1.1 in injA, like 15 partilces

            tmp = convX_mass_xpos(xval,zeroData1[j][k],wp,wUltiVol['C'],wUlit['C'])

            massC_tmp  += tmp[0] * mweightx
            massC0_tmp += tmp[1] * mweightx


    #===================================#
    convX_dict['x'].append(xval)
    #convX_dict['mp'].append(mass_tmp)
    #convX_dict['mChar0'].append(massChar0_tmp)
    convX_dict['mC'].append(massC_tmp)
    convX_dict['mC0'].append(massC0_tmp)
    #convX_dict['mp0'].append(massp0_tmp)

    #if massC0_tmp == 0:
    #    XC = 0
    #else:
        #Xp = (massp0_tmp - mass_tmp)/massp0_tmp
    XC = (massC0_tmp - massC_tmp)/massC0_tmp #* 0.88888

    convX_dict['XC'].append(XC)
    #===================================#

    #print(xval,mass_tmp,massChar0_tmp,massC0_tmp,massp0_tmp,Xp,counter_tmp)
    print(xval,massC_tmp,massC0_tmp,XC)

    #mass_tmp      = 0
    #massChar0_tmp = 0
    massC_tmp     = 0
    massC0_tmp    = 0
    #massp0_tmp    = 0
    #counter_tmp   = 0




print('-- finished --')



