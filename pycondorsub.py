import numpy as np
import os

#DIR="/net/malachit/users/meinecke/diss_cdata/three_section/"
#DIR="/net/granat/users/meinecke/diss_cdata/three_section/"
DIR="/net/granat/users2/meinecke/diss_cdata/three_section/"
PROG="QDMLL_3S"
RUN="44_Jitter_along_U_normal_Noise"

if(os.path.exists(PROG)):
  print "program " + PROG + " found"
  
  PATH = DIR+"run"+RUN
  
  if( os.path.exists(PATH)  ):
    print "Error: directory " + PATH + "already exists"
  
  else:
     
    os.mkdir(PATH)
    os.mkdir(PATH+"/data")
    os.mkdir(PATH+"/data/TS")
    os.mkdir(PATH+"/data/bin")
    os.mkdir(PATH+"/data/powerSpec")
    os.mkdir(PATH+"/data/IntAC")
    os.mkdir(PATH+"/data/posTAvDyn")
    os.mkdir(PATH+"/data/LTTJ")
    os.mkdir(PATH+"/data/LTTJ/powerspecs")
    
    print "Directory " + PATH + " has been created"
    
    argfile = open('pyargfile','w')
    
    
    ###Long long TS for high res spectra
    #MEM=1
    #for r in range(1,1001):
      #argfile.write("-HRPowerSpec -intTime 0.752e6 -outTime -1 -SampleDt 0.25 -maxOutFreq 700")
      #argfile.write(" -noiseStr 1E-6")
      #argfile.write(" -HRS_buffer")
      ##argfile.write(" -HRPowerSpec_Hann")
      #argfile.write(" -strSuffix _run_" +str(r))
      #argfile.write(" -loadHist -histFile " + DIR + "Xhist_stdFML.bin" )
      #argfile.write("\n")
    
    ##LTTJ Long wPS along noiseStr
    #MEM=0
    #for nStr in np.geomspace(1E-17,1E-4,14):
      #for b in np.arange(51,101):
        #argfile.write("-LTTJ -LTTJ_buffer -LTTJ_BIntTime 1e5 -LTTJ_intTime 5e6 -LTTJ_nNoiseRel 1")
        #argfile.write(" -noiseStr " + "{:.6E}".format(nStr) )
        ##argfile.write(" -LTTJ_wPS -LTTJ_PSOut 280 -SampleDt 0.25")
        #argfile.write(" -P_G 0.7 -U 4.0")
        #argfile.write(" -LTTJ_str_suf _nStr_" + "{:.6E}".format(nStr) +"_batch_" + str(b) )
        #argfile.write(" -loadHist -histFile " + DIR + "Xhist_stdFML.bin" )
        #argfile.write("\n")
  
    ##LTTJ along P_G and U
    #MEM=1
    #for PG in np.arange(0.55,0.97,0.00525):
      #for U in np.arange(2,10.001,0.1):
        #argfile.write("-LTTJ -LTTJ_buffer -LTTJ_BIntTime 1e5 -LTTJ_intTime 8e5 -LTTJ_nNoiseRel 25")
        #argfile.write(" -P_G " + "{:1.2f}".format(PG) )
        #argfile.write(" -U " + "{:1.3f}".format(U) )
        #argfile.write(" -noiseStr 1E-14")
        #argfile.write(" -LTTJ_testFML")
        #argfile.write(" -LTTJ_str_suf _P_G_" + "{:1.2f}".format(PG) + "_U_" + "{:1.2f}".format(U) )
        #argfile.write(" -loadHist -histFile " + DIR + "Xhist_stdFML.bin" )
        #argfile.write("\n")
  
  
    ##LTTJ Long wPS along P_G
    #MEM=2
    ##for PG in np.arange(0.78,1.051,0.01):
    #for PG in np.arange(0.78,0.911,0.01):
      #for b in np.arange(1,101):
        #argfile.write("-LTTJ LTTJ_wPS -LTTJ_PSOut 280 -LTTJ_buffer -LTTJ_BIntTime 1e5 -LTTJ_intTime 5e6 -LTTJ_nNoiseRel 1")
        #argfile.write(" -P_G " + "{:1.2f}".format(PG) )
        #argfile.write(" -noiseStr 1E-14")
        #argfile.write(" -SampleDt 0.25")
        #argfile.write(" -LTTJ_str_suf _P_G_" + "{:1.2f}".format(PG) +"_batch_" + str(b) )
        #argfile.write(" -loadHist -histFile " + DIR + "Xhist_stdFML.bin" )
        #argfile.write("\n")
        
        
    ##LTTJ Long wPS along P_G
    #MEM=2
    #for U in np.linspace(5.0,9.0,13.0):
      #for b in np.arange(1,101):
        #argfile.write("-LTTJ LTTJ_wPS -LTTJ_PSOut 280 -LTTJ_buffer -LTTJ_BIntTime 1e5 -LTTJ_intTime 5e6 -LTTJ_nNoiseRel 1")
        #argfile.write(" -U " + "{:1.2f}".format(U) )
        ##argfile.write(" -noiseStr 1E-14")
        #argfile.write(" -SampleDt 0.25")
        #argfile.write(" -LTTJ_str_suf _U_" + "{:1.2f}".format(U) +"_batch_" + str(b) )
        #argfile.write(" -loadHist -histFile " + DIR + "Xhist_stdFML.bin" )
        #argfile.write("\n")
    
        
    ##LTTJ Long wPS along SAStart
    #MEM=2
    #for SAS in np.arange(0.0,0.071,0.005):
      #for b in np.arange(1,101):
        #argfile.write("-LTTJ -LTTJ_wPS -LTTJ_PSOut 280 -LTTJ_buffer -LTTJ_BIntTime 1e5 -LTTJ_intTime 5e6 -LTTJ_nNoiseRel 1")
        #argfile.write(" -SAStart " + "{:1.3f}".format(SAS) )
        #argfile.write(" -noiseStr 1E-6")
        #argfile.write(" -SampleDt 0.25")
        #argfile.write(" -LTTJ_str_suf _SAStart_" + "{:1.3f}".format(SAS) +"_batch_" + str(b) )
        ##argfile.write(" -loadHist -histFile " + DIR + "Xhist_stdFML.bin" )
        #argfile.write("\n")
        
    ##LTTJ Long wPS along theta and PG
    #PVals = np.linspace(0.5, 1.05, 12)
    #ThetaVals = ((PVals-0.53)*4.14 + 0.82)/2.0
    #MEM=2
    #for k in range(PVals.size):
      #for b in np.arange(1,101):
        #argfile.write("-LTTJ -LTTJ_wPS -LTTJ_PSOut 280 -LTTJ_buffer -LTTJ_BIntTime 5e4 -LTTJ_intTime 5e6 -LTTJ_nNoiseRel 1")
        #argfile.write(" -P_G " + "{:1.6f}".format(PVals[k]) )
        #argfile.write(" -theta " + "{:1.6f}".format(ThetaVals[k]) )
        #argfile.write(" -noiseStr 1E-6")
        #argfile.write(" -SampleDt 0.25")
        #argfile.write(" -LTTJ_str_suf _P_G_" + "{:1.6f}".format(PVals[k]) +"_theta_" + "{:1.6f}".format(ThetaVals[k]) + "_batch_" + str(b) )
        #argfile.write(" -loadHist -histFile " + DIR + "Xhist_stdFML.bin" )
        #argfile.write("\n")
  
    ##LTTJ long
    #MEM=4
    #for b in np.arange(1,251):
      #argfile.write("-LTTJ -LTTJ_buffer -LTTJ_BIntTime 5e4 -LTTJ_intTime 1e7 -LTTJ_nNoiseRel 4")
      #argfile.write(" -LTTJ_str_suf _batch_" + str(b) )
      #argfile.write(" -LTTJ_wPS -LTTJ_PSOut 700" )
      #argfile.write(" -noiseStr 1E-6")
      #argfile.write(" -loadHist -histFile " + DIR + "Xhist_stdFML.bin" )
      #argfile.write("\n")
      
      
    ##U
    #MEM=1
    #for U in np.linspace(0,10,201, endpoint=True):
      #argfile.write("-sweep_P_G_wU -U " + "{:1.3f}".format(U) +" -sIntTime 120000 -sOutTime 40000 -sStart 0.4 -sEnd 1.6 -sSteps 200")
      ##argfile.write(" wpowerSpec wIntAC wExtrema wDeltaT")
      #argfile.write(" -noiseStr 1E-6")
      #argfile.write(" -wNoSweep")
      ##argfile.write(" -loadHist -histFile " + DIR + "Xhist_stdHML3.bin")
      ##argfile.write(" -loadHist -histFile " + DIR + "Xhist_stdFML.bin")
      #argfile.write(" -loadHist -histFile " + DIR + "Xhist_P12_GESQ4.bin")
      #argfile.write("\n")
      
    ##U for mean
    #MEM=1
    #for b in range(1,11):
      #argfile.write("-sweep_P_G_wU -U 6.0 -sIntTime 100000 -sOutTime 50000 -sStart 0.4 -sEnd 1.6 -sSteps 400")
      #argfile.write(" -strSuffix _run_"+str(b))
      #argfile.write(" -noiseStr 1E-6")
      #argfile.write(" wpowerSpec -wIntAC wExtrema wNoSweep wDeltaT")
      ##argfile.write(" -loadHist -histFile " + DIR + "Xhist_stdHML3.bin")
      #argfile.write("\n")
      
    ##av specs mean
    #MEM=1
    #for b in range(1,21):
      #argfile.write("-sweep_P_G_wU -U 4.0 -sIntTime 90000 -sOutTime 30000 -sStart 1.5 -sEnd 1.5 -sSteps 1")
      #argfile.write(" -strSuffix _run_"+str(b))
      #argfile.write(" -noiseStr 1E-6")
      #argfile.write(" -wpowerSpec wIntAC wExtrema wNoSweep wDeltaT")
      #argfile.write(" -loadHist -histFile " + DIR + "Xhist_U4.bin")
      #argfile.write("\n")
      
    ##P_G with noiseStr
    #MEM=1
    #for nStr in np.geomspace(1E-17,1E-4,201):
      #argfile.write("-sweep_P_G_wNoiseStr -noiseStr " + "{:.6E}".format(nStr) + " -sIntTime 120000 -sOutTime 60000 -sStart 0.4 -sEnd 1.6 -sSteps 200")
      ##argfile.write(" wpowerSpec wIntAC wExtrema wDeltaT")
      #argfile.write(" -wNoSweep")
      ##argfile.write(" -loadHist -histFile " + DIR + "Xhist_stdHML3.bin")
      #argfile.write("\n")
      
      
    ##SAStart
    #MEM=1
    ##for SAS in np.linspace(0,0.07,201, endpoint=True):
    #for SAS in np.linspace(0,0.040250,116, endpoint=True):
      #argfile.write("-sweep_P_G_wSAStart -SAStart " + "{:1.6f}".format(SAS) +" -sIntTime 120000 -sOutTime 40000 -sStart 0.6 -sEnd 1.4 -sSteps 200")
      #argfile.write(" -wNoSweep")
      ##argfile.write(" -loadHist -histFile " + DIR + "Xhist_stdFML.bin")
      ##argfile.write(" -loadHist -histFile " + DIR + "Xhist_stdHML3.bin")
      ##argfile.write(" -loadHist -histFile " + DIR + "Xhist_P12_GESQ4.bin")
      #argfile.write(" -loadHist -histFile " + DIR + "Xhist_FML_081_SAS_002.bin")
      ##argfile.write(" wpowerSpec wIntAC wExtrema wNoSweep -wDeltaT")
      #argfile.write("\n")
      
      
      
    ##SAStart for quasi const PG
    #MEM=1
    #for SAS in np.linspace(0.05,0.07,201, endpoint=True):
      #argfile.write("-sweep_P_G_wSAStart -SAStart " + "{:1.6f}".format(SAS) +" -sIntTime 120000 -sOutTime 60000 -sStart 1.200 -sEnd 1.200100 -sSteps 100")
      #argfile.write(" wpowerSpec -wIntAC -wExtrema -wNoSweep -wDeltaT")
      ##argfile.write(" -loadHist -histFile " + DIR + "Xhist_stdHML3.bin")
      #argfile.write("\n")
      
    ##GESQ
    #MEM=1
    #for gamma in np.geomspace(0.05,7.5,201):
      #argfile.write("-sweep_P_G_wgamma_ES_Q -gamma_ES_Q " + "{:.6E}".format(gamma) +" -wExtrema -wNoSweep -sIntTime 120000 -sOutTime 40000 -sStart 0.4 -sEnd 1.6 -sSteps 200")
      #argfile.write(" -noiseStr 1E-6")
      #argfile.write(" -wNoSweep")
      ##argfile.write(" -loadHist -histFile " + DIR + "Xhist_stdHML3.bin")
      #argfile.write(" -loadHist -histFile " + DIR + "Xhist_P12_GESQ4.bin")
      #argfile.write("\n")
      
    ##SARS
    #MEM=1
    #for SARS in np.linspace(0.0,4.0,201, endpoint=True):
      #argfile.write("-sweep_P_G_wSARS -SARS " + "{:1.6f}".format(SARS) +" -sIntTime 120000 -sOutTime 40000 -sStart 0.4 -sEnd 1.6 -sSteps 200")
      #argfile.write(" -noiseStr 1E-6")
      #argfile.write(" -wNoSweep")
      ##argfile.write(" -loadHist -histFile " + DIR + "Xhist_stdFML.bin")
      ##argfile.write(" -loadHist -histFile " + DIR + "Xhist_stdHML3.bin")
      #argfile.write(" -loadHist -histFile " + DIR + "Xhist_P12_GESQ4.bin")
      #argfile.write("\n")
      
    ##Theta
    #MEM=1
    #for SAS in np.linspace(0.0,1.5,201, endpoint=True):
      #argfile.write("-sweep_P_G_wtheta -theta " + "{:1.6f}".format(SAS) +" -sIntTime 120000 -sOutTime 40000 -sStart 0.4 -sEnd 1.6 -sSteps 200")
      #argfile.write(" -noiseStr 1E-6")
      #argfile.write(" -wNoSweep")
      ##argfile.write(" -loadHist -histFile " + DIR + "Xhist_stdFML.bin")
      ##argfile.write(" -loadHist -histFile " + DIR + "Xhist_stdHML3.bin")
      ##argfile.write(" -loadHist -histFile " + DIR + "Xhist_FML_theta_041_P0525.bin")
      ##argfile.write(" -loadHist -histFile " + DIR + "Xhist_HML3_theta015_P055.bin")
      #argfile.write(" -loadHist -histFile " + DIR + "Xhist_P12_GESQ4.bin")
      #argfile.write("\n")
      
    ##Theta zoom
    #MEM=1
    #for SAS in np.linspace(0.1,0.5,161, endpoint=True):
      #argfile.write("-sweep_P_G_wtheta -theta " + "{:1.6f}".format(SAS) +" -sIntTime 80000 -sOutTime 30000 -sStart 0.425 -sEnd 0.575 -sSteps 150")
      #argfile.write(" -wExtrema -wNoSweep")
      #argfile.write(" -loadHist -histFile " + DIR + "Xhist_stdFML_theta0.5.bin")
      #argfile.write("\n")
      
      
    ##LTTJ with -posTAvDyn
    #MEM=14
    #for b in range(2000):
      #argfile.write( "-posTAvDyn -intTime 12e4 -outTime 8e4")
      #argfile.write(" -noiseStr 1E-6")
      #argfile.write(" -SAStart 0.01")
      #argfile.write(" -strSuffix _run_"+str(b))
      #argfile.write(" -loadHist -histFile " + DIR + "Xhist_stdSAleftFML.bin")
      #argfile.write("\n")
      
    ##-TS for mean optical spec
    #MEM=4
    #for b in range(100):
      #argfile.write( "-TS -intTime 15e4 -outTime 10e4")
      #argfile.write(" -noiseStr 1E-6")
      #argfile.write(" -strSuffix _run_"+str(b))
      #argfile.write(" -opticalSpec")
      #argfile.write(" -P_G 1.08")
      #argfile.write(" -loadHist -histFile " + DIR + "Xhist_stdHML3.bin")
      #argfile.write("\n")
      
    
    ##posTAvDyn along Theta
    #MEM=16
    #PVals = np.linspace(0.5, 1.05, 201)
    #ThetaVals = ((PVals-0.53)*4.14 + 0.82)/2.0
    #for k in range(PVals.size):
      #argfile.write( "-posTAvDyn -intTime 12e4 -outTime 4e4")
      #argfile.write(" -noiseStr 1E-6")
      #argfile.write(" -P_G " + "{:1.6f}".format(PVals[k]) )
      #argfile.write(" -theta " + "{:1.6f}".format(ThetaVals[k]) )
      #argfile.write(" -strSuffix _P_G_"+"{:1.6f}".format(PVals[k]) + "_theta_"+"{:1.6f}".format(ThetaVals[k]) )
      #argfile.write(" -loadHist -histFile " + DIR + "Xhist_stdFML120.bin")
      #argfile.write("\n")
      
    ##posTAvDyn along P
    #MEM=16
    #PVals = np.linspace(0.6, 1.00, 201)
    #for k in range(PVals.size):
      #argfile.write( "-posTAvDyn -intTime 12e4 -outTime 4e4")
      #argfile.write(" -noiseStr 1E-6")
      #argfile.write(" -P_G " + "{:1.6f}".format(PVals[k]) )
      #argfile.write(" -strSuffix _P_G_"+"{:1.6f}".format(PVals[k]) )
      #argfile.write(" -loadHist -histFile " + DIR + "Xhist_stdFML120.bin")
      #argfile.write("\n")
      
      
    ##posTAvDyn along U
    #MEM=16
    #UVals = np.linspace(4.0, 10.0, 201)
    #for k in range(UVals.size):
      #argfile.write( "-posTAvDyn -intTime 12e4 -outTime 4e4")
      #argfile.write(" -noiseStr 1E-6")
      #argfile.write(" -U " + "{:1.6f}".format(UVals[k]) )
      #argfile.write(" -strSuffix _U_"+"{:1.6f}".format(UVals[k]) )
      #argfile.write(" -loadHist -histFile " + DIR + "Xhist_stdFML120.bin")
      #argfile.write("\n")
      
      
    ##posTAvDyn along SAStart
    #MEM=16
    #SASVals = np.linspace(0,0.07,201, endpoint=True)
    #for k in range(SASVals.size):
      #argfile.write( "-posTAvDyn -intTime 12e4 -outTime 4e4")
      #argfile.write(" -noiseStr 1E-6")
      #argfile.write(" -SAStart " + "{:1.6f}".format(SASVals[k]) )
      #argfile.write(" -strSuffix _SAStart_"+"{:1.6f}".format(SASVals[k]) )
      #argfile.write(" -loadHist -histFile " + DIR + "Xhist_stdFML120.bin")
      #argfile.write("\n")
  
    argfile.close()
    print "pyargfile created"
    
    substr = "qsub -mem " + str(MEM) + " -m n -speed 4 -w " + PATH + " -argfile pyargfile " + PROG
    print "Submit string: "
    print substr
    
    os.system( substr )

else:
  print "Error: program not found!"

