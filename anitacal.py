## **************************************************************************
## AnitaEventCal
##  This is a very hacky version of the AnitaEventCalibrator code from https://github.com/anitaNeutrino/eventReaderRoot/tree/master
## **************************************************************************


import numpy as np

class AnitaEventCal:

    def __init__(self):
        self.epsilons=np.zeros((12,4,2))
        self.binValues=np.zeros((12,4,2,260))
        self.timeValues=np.zeros((12,4,258*3))   ##774 samples RCO 0 -> RCO 1 -> RCO 0
        self.mvCalib=np.ones((12,9,4))  #Don't have voltage calibration for clock channels so leave at 1
        self.relCableDelays=np.zeros((12,9,4))
        self.relPhaseCentreDelays=np.zeros((12,9))
        self.readEpsilons()
        self.readBinByBin()
        self.readVoltageCalib()
        self.readRelativeCableDelays()
        self.readRelativePhaseCentreDelays()
        self.fillTimeValues()
        
    
    
    def readVoltageCalib(self):
        f = open("calib/simpleVoltageCalibrationAnita4.txt", "r")
        for line in f.readlines()[1:]:
            tokens=line.split()
            surf=int(tokens[0])
            chan=int(tokens[1])
            chip=int(tokens[2])
            mvperadc=float(tokens[3])
            self.mvCalib[surf][chan][chip]=mvperadc

    
    def readRelativeCableDelays(self):
        f = open("calib/relativeCableDelaysAnita4.dat", "r")
        for line in f.readlines()[1:]:
            tokens=line.split()
            surf=int(tokens[0])
            chan=int(tokens[1])
            chip=int(tokens[2])
            delay=float(tokens[3])
            self.relCableDelays[surf][chan][chip]=delay

    
    def readRelativePhaseCentreDelays(self):
        f = open("calib/relativePhaseCenterToAmpaDelaysAnita4.dat", "r")
        for line in f.readlines()[1:]:
            tokens=line.split()
            surf=int(tokens[0])
            chan=int(tokens[1])
            delay=float(tokens[2])
            self.relPhaseCentreDelays[surf][chan]=delay

    
    def readEpsilons(self):
        f = open("calib/epsilonFromBenS.dat", "r")
        for line in f.readlines()[1:]:
            tokens=line.split()
            surf=int(tokens[0])
            chip=int(tokens[1])
            rco=int(tokens[2])
            binwidth=float(tokens[3])
            self.epsilons[surf][chip][rco]=binwidth

    
    def readBinByBin(self):
        f = open("calib/justBinByBin.dat", "r")
        for line in f.readlines()[1:]:
            tokens=line.split()
            surf=int(tokens[0])
            chip=int(tokens[1])
            rco=int(tokens[2])
            samp=int(tokens[3])
            binwidth=float(tokens[4])
            self.binValues[surf][chip][rco][samp]=binwidth

    
    def fillTimeValues(self):
        # Here we are going to do something with the timing calibration numbers
        # Reminder there are 260 capacitors per LAB chip and the speed the switching happens is depend on the RCO phase (high or low of clock)
        # For now we will throw away the first (0) and last capacitor (259) so we have just 258 capacitors


        for surf in range(12):
            for chip in range(4):
                
                # First step is get all the bin values for RCO=0 and form the cumulative sum
                bv0=np.cumsum(np.concatenate ( ([0],self.binValues[surf][chip][0][:-1]) ))
                # Next step is get all the bin values for RCO=1 and form the cumulative sum
                bv1=np.cumsum(np.concatenate ( ([0],self.binValues[surf][chip][1][:-1]) ))
                # Now we can form basically a ring buffer of times RCO0 -> RCO1 -> RCO0
                # Wherever our event actually lies we can then pick the appropriate times and zero adjust
                bvBig=np.concatenate( (bv0[1:-1],bv1[1:-1]+bv0[255]+self.epsilons[surf][chip][1]))
                bvBig=np.concatenate( (bvBig, bvBig[-4]+bv0[1:-1]+self.epsilons[surf][chip][0]))
                
                #print(np.shape(bvBig))
                #print(np.diff(bvBig))
                
                #Should add this as a check somewhere in the calibration code
                if not np.all(np.diff(bvBig) > 0):
                    raise ValueError("Some of our time differences are negative.... this shouldn't happen")
                self.timeValues[surf][chip]=np.array(bvBig)


    def getLabChip(self,chipIdFlag, chanIndex):
        return chipIdFlag[chanIndex]&0x3 # Returns the LABRADOR number
    def getRCO(self,chipIdFlag, chanIndex):
        return (chipIdFlag[chanIndex]&0x4)>>2 #Returns the RCO phase
    def getFirstHitBus(self,firstHitbus, chanIndex):
        return firstHitbus[chanIndex] #Returns the firstHitbus value for the channel
    
    def getLastHitBus(self,lastHitbus, chipIdFlag, chanIndex): # Returns the lastHitbus value for the channel
        if(lastHitbus[chanIndex]<255): 
            return lastHitbus[chanIndex]
        return (lastHitbus[chanIndex]) + ((chipIdFlag[chanIndex])>>4) 

    def  getWrappedHitBus(self,chipIdFlag, chanIndex) : #  Return the wrapped hitbus flag for the channel. When the HITBUS is wrapped the waveform runs from firstHitbus+1 to lastHitbus-1, otherwise it runs from lastHitbus+1 to firstHitbus-1 (crossing the 259-->0 boudnary).
        return ((chipIdFlag[chanIndex])&0x8)>>3

    def getEarliestSample(self,fh,lh,wh):
        #lh=getLastHitBus(lastHitbus,chipIdFlag,chanIndex)
        #fh=getFirstHitBus(firstHitbus,chanIndex)
        #wh=getWrappedHitBus(chipIdFlag,chanIndex)
        earliestSample=0
        if not wh:
            earliestSample=lh+1
        else:    
            earliestSample=fh+1
     
        if(earliestSample==0):
            earliestSample=1
        if(earliestSample<259):
            return earliestSample;
        return 1;
    
    
    def getLatestSample(self,fh,lh,wh):
        #lh=getLastHitBus(lastHitbus,chipIdFlag,chanIndex)
        #fh=getFirstHitBus(firstHitbus,chanIndex)
        #wh=getWrappedHitBus(chipIdFlag,chanIndex)
        latestSample=258
        if not wh: 
            latestSample=fh-1
        else:
           latestSample=lh-1
         
        if(latestSample>0): 
            return latestSample;
        return 258;


    def getIndicesForSurfs(self,eventDict,calDict):
        validTimeInds=np.zeros((12,2),dtype=int)
        
        rco=calDict['rcoArray']
        cifList=eventDict['chipIdFlag[108]']
        fhList=eventDict['firstHitbus[108]']
        lhList=eventDict['lastHitbus[108]']
        for surf in range(12):
            for chan in range(1):
                chanIndex=9*surf + chan
                fh=self.getFirstHitBus(fhList,chanIndex)
                lh=self.getLastHitBus(lhList,cifList,chanIndex)
                wh=self.getWrappedHitBus(cifList,chanIndex)
                earliestSample=self.getEarliestSample(fh,lh,wh)
                latestSample=self.getLatestSample(fh,lh,wh)
                endRco=np.copy(rco[surf])
                startRco=np.copy(rco[surf])
                if(earliestSample>latestSample):
                    startRco=1-startRco
                #print(surf,startRco,endRco)
                    
                startIndex=int((258*startRco)+earliestSample)
                endIndex=int((258*startRco)+(258*abs(startRco-endRco))+latestSample)
                #print("fhList",fhList)
                #print("lhList",lhList)
                #print("fh",fh)
                #print("lh",lh)
                #print(fh,lh,wh,earliestSample,latestSample,startIndex,endIndex)
                
                if(endIndex>startIndex+250):
                    endIndex=startIndex+250

                #if(endIndex<startIndex+240):
                #    raise ValueError("For some reason we only have "+str(endIndex-startIndex)+" samples for SURF "+str(surf))  
                validTimeInds[surf][0]=startIndex
                validTimeInds[surf][1]=endIndex
        return validTimeInds
    
    
    def calibrateEvent(self,eventDict,calDict):

        #//! The plan for calibration happiness:
        #//! Step 1: figure out from WaveCalType_t calType exactly what we're going to do  -- RJNPyHack Will assume full calibration
        #//! Step 2: Remove the spiky clock
        #//! Step 3: Figure out RCO phase from clock period (if bin-to-bin dts required)
        #//! Step 4: Update rolling temperature correction (+copy to event)
        #//! Step 5: Unwrap all channels (if requested)
        #//! Step 6: Apply bin-to-bin timing (if requested)
        #//! Step 7: Apply voltage calib (if requested)
        #//! Step 9: Zero mean all non-clock channels
        ##//! Step 8: Find trigger jitter correction
        #//! Step 10: Apply channel-to-channel cable delays
        #//! Step 11: Copy everything to UsefulAnitaEvent

        # Initial states for calibration flags, modified based on calType below
        fNominal                        = False;
        fApplyTempCorrection            = True; 
        fRemoveClockSpike = True;
        fUnwrap = True;
        fBinToBinDts = True;
        fVoltage = True;
        fApplyTriggerJitterCorrection = True;
        fApplyCableDelays = True;
        fZeroMeanNonClockChannels = True;
        fApplyExtraDelayFromPhaseCenter = True;

        validTimeInds=self.getIndicesForSurfs(eventDict,calDict)
        #print(validTimeInds)

        adcVals=np.zeros((12,9,250),dtype=int) #In the above calibration scheme 250 are the maximumn number of valid samples
        for surf in range(12):
            si=validTimeInds[surf][0]%258
            ei=validTimeInds[surf][1]%258
            for chan in range(9):
                chanIndex=9*surf+chan
                if ei>si:
                    adcVals[surf][chan][0:ei-si]=eventDict['data[108][260]'][chanIndex][si+1:ei+1]
                else:
                    adcVals[surf][chan][0:258-si]=eventDict['data[108][260]'][chanIndex][si+1:259]
                    adcVals[surf][chan][258-si:(258-si)+ei]=eventDict['data[108][260]'][chanIndex][1:ei+1]

        
        #//! Step 9: Zero mean all non-clock channels
        adcOffset=np.zeros((12,9))
        for surf in range(12):
            for chan in range(8): #Excluding clock channels
                adcOffset[surf][chan]=np.mean(adcVals[surf][chan])

        chips=np.zeros((12),dtype=int)
        timeScale=np.ones((12,9))
        timeOffset=np.zeros((12,9))
        mvScale=np.ones((12,9))  #Okay so this is fixed and does not vary event by event
        for surf in range(12):
            dt=-1*calDict['clockPhiArray'][surf]
            tscale=calDict['tempFactorGuesses'][surf]
            chanIndex=9*surf
            chips[surf]=self.getLabChip(eventDict['chipIdFlag[108]'],chanIndex)
            for chan in range(9):
                chanIndex=9*surf+chan
                chip=self.getLabChip(eventDict['chipIdFlag[108]'],chanIndex)
                mvScale[surf][chan]=self.mvCalib[surf][chan][chip]
                timeScale[surf][chan]=tscale
                timeOffset[surf][chan]=dt+self.relCableDelays[surf][chan][chip]-self.relPhaseCentreDelays[surf][chan]


        
        return adcVals,adcOffset,mvScale,timeOffset,timeScale,validTimeInds,chips


        #Step 2: For now skip the spiky clock
        #Step 3: We will have fRcoArray and fClockProblem from the Calibrated Info files
        #Step 4: We will fave fTempFactorGuesses from Calibrated Info files
        #Step 5: Now we do the unwrapping


        
