# BEGIN PLOT /MC_ColourReconstruction/d*
RatioPlotYMin=0.6
RatioPlotYMax=1.4
RatioPlotYLabel=Ratio
XTwosidedTicks=1
YTwosidedTicks=1
LogY=0
LeftMargin=1.8
RightMargin=4.5
Variations=fsr_murfac=0.5,fsr_murfac=2.0,isr_murfac=0.5,isr_murfac=2.0,
BandComponentEnv=fsr_murfac=0.5,fsr_murfac=2.0,isr_murfac=0.5,isr_murfac=2.0,
# END PLOT

# BEGIN PLOT /MC_ColourReconstruction/ccjetsinvariantmass
Title=Mass of $W$ boson calculated from jets identified as colour connected 
XLabel=$M_{W}$ GeV
YLabel=$\dfrac{1}{\sigma}\dfrac{\text{d}\sigma}{\text{d}M_W}$
# END PLOT

# BEGIN PLOT /MC_ColourReconstruction/lightjetsinvariantmass
Title=Mass of $W$ boson calculated from all light jet combinations
XLabel=$M_{W}$ GeV
YLabel=$\dfrac{1}{\sigma}\dfrac{\text{d}\sigma}{\text{d}M_W}$
# END PLOT

# BEGIN PLOT /MC_ColourReconstruction/leading_lightjetsinvariantmass
Title=Mass of $W$ boson calculated from leading light jets 
XLabel=$M_{W}$ GeV
YLabel=$\dfrac{1}{\sigma}\dfrac{\text{d}\sigma}{\text{d}M_W}$
# END PLOT

# BEGIN PLOT /MC_ColourReconstruction/numberwjetpairs
Title=Number of light jet pair combinations that have invariant mass within $40$ GeV of $W$ mass 
XLabel=Number of pairs
YLabel=Count
# END PLOT

# BEGIN PLOT /MC_ColourReconstruction/numberccjetpairs
Title=Number of light jet pair combinations that are identified as colour connected
XLabel=Number of pairs
YLabel=Count
# END PLOT

# BEGIN PLOT /MC_ColourReconstruction/pull_12
Title=(forwards) Pull angle between two light jets 
YMin=0.8
XLabel=$\phi(j_1, j_2) / \pi$
YLabel=$\dfrac{1}{\sigma}\dfrac{\text{d}\sigma}{\text{d}\phi(j_1,j_2)}$
# END PLOT
