##import numpy as np

origFNameL = {'lecco':'lecco'}
fNameStep = {'lecco': 0.1667}

#origFNameL = ['carenno']
defVal = '0.0'
addIntDays = 2
timeIntMin = 10
minHour = 60
hourDay = 24
minDay = minHour*hourDay
stepDay = int(minDay/timeIntMin)
iniDay = 26 # originally in ev 26
iniMnt = 10 # originally in ev 10
iniYer = 2018 # originally in ev 2018
### Search for a date starting
print("stepDay", stepDay,addIntDays)
for oFNL in origFNameL:
  step = fNameStep[oFNL]
  step = int(step *3600)
  oFNLv = origFNameL[oFNL] #variation for the internal name
  nameF = (oFNL+'/temperature_'+oFNLv+'_ev.txt')
  outF  = (oFNL+'/temperature_'+oFNLv+'_ev0.txt')
  fOut = open(outF,'w')
  hourL = 3600
  with open(nameF,'r') as f:
    ll = f.readline() ## Header line
    fOut.write(ll)
    ll = f.readline() # read second line
    lsep = ll.strip().split()
    idFlag = lsep[0]
    preT = float(lsep[3])
    preTC = "{:2.1f}".format(preT)
    print(preTC)
    write_ev = False
    for ii in range(22,iniDay):
      #print("puttanamadonna", iniMnt,iniYer)
      date = "{0:02d}".format(int(ii))
      date1 = "{0:02d}".format(int(iniMnt))
      date2 = "{0:04d}".format(int(iniYer))
      dateTot =date2 +"/"+date1+"/"+date
      #print(date,date1,date2,dateTot)
      ## time
      #print("max n step day", oFNL, int(24*hourL/step))
      for jj in range(0,int(24*hourL/step)):
        timeH = int(step*jj)/hourL
        timeM = (step*jj)%hourL/60
        #print(timeM,timeH)
        timeStr0 = "{0:02d}".format(int(timeH))
        timeStr1 = "{0:02d}".format(int(timeM))
        timeStr = timeStr0 + ":"+ timeStr1
        #print(date, timeStr,preRain)
        ## Id Day Hour Value
        fOut.write("{} {} {} {}\n".format(idFlag,dateTot,timeStr,preTC))
      #break
    while ll :
      lsep = ll.strip().split()
      #print(lsep)
      idFlag = lsep[0]
      currDate = lsep[1].split("/")
      day = int(currDate[2])
      #if(int(currDate[0]) == iniYer):
      #  if(int(currDate[1]) == iniMnt):
      #    if(int(currDate[2]) == iniDay):
      write_ev = True
      if write_ev:
        fOut.write(ll)
      ll = f.readline() # read second line
      

    #iniDay = day - addIntDays
    ## Filling the additional lines
    #for ii in range(0,addIntDays*stepDay):
    #for ii in range(0,13):
    #  totTime = ii * timeIntMin
    #  addDay  = int(totTime/minDay)
    #  addHour = int((totTime%minDay)/minHour)
    #  addMin  = (totTime%minDay)%minHour
    #  hSign = "{:02d}:{:02d}".format(addHour,addMin)
    #  dSign = "{}/{}/{:02d}".format(currDate[0],currDate[1],day-addIntDays+addDay)
      
    #  #print(totTime, addDay,addHour,addMin, hSign,dSign)
    #  print("{} {} {} {}".format(idFlag,dSign,hSign,defVal))
    
