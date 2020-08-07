# BEGIN PLOT /MC_ColourFlow/d*
RatioPlotYMin=0.6
RatioPlotYMax=1.4
RatioPlotYLabel=Ratio
XTwosidedTicks=1
YTwosidedTicks=1
LogY=0
LeftMargin=1.8
RightMargin=3
Variations=fsr_murfac=0.5,fsr_murfac=2.0
BandComponentEnv=fsr_murfac=0.5,fsr_murfac=2.0
# END PLOT

# BEGIN PLOT /MC_ColourFlow/pull_12_all
Title=Pull angle 1,2 between $W$-jets, all constituents
RatioPlotYMin=0.99
RatioPlotYMax=1.01
YMin=0.8
XLabel=$\phi(j^W_1, j^W_2) / \pi$
YLabel=$\dfrac{1}{\sigma}\dfrac{\text{d}\sigma}{\text{d}\phi(j^W_1,j^W_2)}$
# END PLOT


# BEGIN PLOT /MC_ColourFlow/pull_21_all
Title=Pull angle 2,1 between $W$-jets, all constituents
RatioPlotYMin=0.95
RatioPlotYMax=1.05
YMin=0.8
XLabel=$\phi(j^W_2, j^W_1) / \pi$
YLabel=$\dfrac{1}{\sigma}\dfrac{\text{d}\sigma}{\text{d}\phi(j^W_1,j^W_2)}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/pull_both_all
Title=Pull angle (1,2 and 2,1) between $W$-jets, all constituents
RatioPlotYMin=0.90
RatioPlotYMax=1.10
YMin=0.8
XLabel=$\phi(j^W_1, j^W_2) / \pi + \phi(j^W_2, j^W_1) / \pi$
YLabel=$\dfrac{1}{\sigma}\dfrac{\text{d}\sigma}{\text{d}\phi(j^w_1,j^w_2)}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/pull_angle_pullmagnitude_cor
Title=Correlation between pull angle pull magnitude
YMin=0.002
XLabel=$\phi(j^W_1, j^W_2) / \pi$
YLabel=$|t|$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/pull_angle_asymm
Title=Asymmetry of pull angle 1,2 between $W$-jets
XLabel=$\phi(j^W_1, j^W_2) / \pi$
YLabel=$A_{\phi}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/perpproj_asymm
Title=Perpendicular asymmetry of pull vector of highest p_t $W$-jet
XLabel=$|t_{\bot}|$
YLabel=$A_{\bot}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/paraproj_asymm
Title=Parallel asymmetry of pull vector of highest p_t $W$-jet
RatioPlotYMin=0.8
RatioPlotYMax=1.2
XLabel=$|t_{\parallel}|$
YLabel=$A_{\parallel}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/paraproj_asymm_onebin
Title=Parallel asymmetry of pull vector of highest p_t $W$-jet
XLabel=$|t_{\parallel}|$
YLabel=$A_{\parallel}$
RatioPlotYMin=0.7
RatioPlotYMax=1.3
# END PLOT

# BEGIN PLOT /MC_ColourFlow/perpproj_asymm_onebin
Title=Perpendicular asymmetry of pull vector of highest p_t $W$-jet
XLabel=$|t_{\bot}|$
YLabel=$A_{\bot}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/perpproj_asymm_M
Title=(weighted) Perpendicular asymmetry of pull vector of highest p_t $W$-jet
XLabel=$|t_{\bot}|$
YLabel=$A_{\bot}^M$ 
# END PLOT

# BEGIN PLOT /MC_ColourFlow/paraproj_asymm_M
Title=(weighted) Parallel asymmetry of pull vector of highest p_t $W$-jet
XLabel=$|t_{\parallel}|$
YLabel=$A_{\parallel}^M$ 
# END PLOT

# BEGIN PLOT /MC_ColourFlow/pull_angle_paraproj_asymm_cor
Title=Correlation between pull angle (mod $\pi$/2) and parallel asymmetry
XLabel=$\phi(j^W_1, j^W_2) / \pi$
YLabel=$A_{\parallel}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/paraproj_pull_angle_cor
Title=Correlation between parallel pull projection and pull angle
XLabel=$t_{\parallel}$
YLabel=$\phi(j^W_1, j^W_2) / \pi$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/paraproj_pull_angle_cor_signed
Title=(signed) Correlation between parallel pull projection and pull angle
XLabel=$t_{\parallel}$
YLabel=$\phi(j^W_1, j^W_2) / \pi$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/paraproj_pullmagnitude_cor
Title=Correlation between parallel pull projection and pull magnitude
XLabel=$|t_{\parallel}|$
YLabel=$|t|$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/pullmagnitude_paraproj_asymm_cor
Title=Correlation between pull magnitude and parallel asymmetry
XLabel=$|t|$
YLabel=$A_{\parallel}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/paraproj_asymm_tail
Title=Parallel asymmetry of pull vector of highest p_t $W$-jet
XLabel=$|t_{\parallel}|$
YLabel=$A_{\parallel}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/bpull_12_all
Title=Pull angle 1,2 between $B$-jets, all constituents
YMin=0.8
YMax=1.1
XLabel=$\phi(j^B_1, j^B_2) / \pi$
YLabel=$\dfrac{1}{\sigma}\dfrac{\text{d}\sigma}{\text{d}\phi(j^B_1,j^B_2)}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/bpull_21_all
Title=Pull angle 2,1 between $B$-jets, all constituents
YMin=0.8
YMax=1.1
XLabel=$\phi(j^B_2, j^B_1) / \pi$
YLabel=$\dfrac{1}{\sigma}\dfrac{\text{d}\sigma}{\text{d}\phi(j^B_2,j^B_1)}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/bparaproj_asymm
Title=Parallel asymmetry of pull vector of highest p_t $B$-jet
XLabel=$|t_{\parallel}|$
YLabel=$A_{\parallel}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/bperpproj_asymm
Title=Perpendicular asymmetry of pull vector of highest p_t $B$-jet
XLabel=$|t_{\bot}|$
YLabel=$A_{\bot}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/pull_angle_12_jetpTbins
Title=Count of pull angle $1,2 < 1/3$ $\pi$ between $W$-jets in jet 1 pT bins
RatioPlotYMin=0.95
RatioPlotYMax=1.05
XLabel=$\pT$ (GeV)
YLabel=$\dfrac{1}{\sigma}\dfrac{\text{d}\sigma}{\text{d}\phi(j^W_1,j^W_2)}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/pull_angle_21_jetpTbins
Title=Count of pull angle $2,1 < 1/3$ $\pi$ between $W$-jets in jet 2 pT bins
RatioPlotYMin=0.95
RatioPlotYMax=1.05
XLabel=$\pT$ (GeV)
YLabel=$\dfrac{1}{\sigma}\dfrac{\text{d}\sigma}{\text{d}\phi(j^W_2,j^W_1)}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/paraproj_asymm_j2
Title=Parallel asymmetry of pull vector $W$-jet 2
XLabel=$|t_{\parallel}|$
YLabel=$A_{\parallel}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/paraproj_asymm_jetpTbins
Title=Parallel asymmetry of pull vector of highest p_t $W$-jet in jet 1 pT bins, $|t_{\parallel}| < 0.0033$
XLabel=$\pT$ (GeV)
YLabel=$A_{\parallel}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/paraproj_asymm_j2_jetpTbins
Title=Parallel asymmetry of pull vector of $W$-jet 2 in jet 2 pT bins, $|t_{\parallel}| < 0.0033$
XLabel=$\pT$ (GeV)
YLabel=$A_{\parallel}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/paraproj_asymm_j2_onebin
Title=Parallel asymmetry of pull vector of $W$-jet 2
XLabel=$|t_{\parallel}|$
YLabel=$A_{\parallel}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/jetpT_pullmagnitude_cor
Title=Correlation between jetpT and $|t|$
XLabel=$\pT$ (GeV)
YLabel=$|t|$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/jetpT_paraproj_cor
Title=Correlation between jetpT and $|t_{\parallel}|$
XLabel=$\pT$ (GeV)
YLabel=$|t_{\parallel}|$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/pull_12_all_lowpT
Title=Pull angle 1,2 between $W$-jets, all constituents, jetpT $<$ 60 GeV
RatioPlotYMin=0.95
RatioPlotYMax=1.05
YMin=0.8
XLabel=$\phi(j^W_1, j^W_2) / \pi$
YLabel=$\dfrac{1}{\sigma}\dfrac{\text{d}\sigma}{\text{d}\phi(j^W_1,j^W_2)}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/paraproj_asymm_lowpT
Title=Parallel asymmetry of pull vector of highest p_t $W$-jet, jetpT $<$ 60 GeV
RatioPlotYMin=0.8
RatioPlotYMax=1.2
XLabel=$|t_{\parallel}|$
YLabel=$A_{\parallel}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/paraproj_asymm_onebin_jetpTbins
Title=Parallel asymmetry of pull vector of highest p_t $W$-jet in jet 1 pT bins
XLabel=$\pT$ (GeV)
YLabel=$A_{\parallel}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/paraproj_asymm_j2_onebin_jetpTbins
Title=Parallel asymmetry of pull vector of $W$-jet 2 in jet 2 pT bins
XLabel=$\pT$ (GeV)
YLabel=$A_{\parallel}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/wJet1_pT
Title=pT of W-jet 1 
XLabel=$\pT$ (GeV)
#YLabel=$\dfrac{1}{\sigma}\dfrac{\text{d}\sigma}{\text{d}\pT$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/wJet2_pT
Title=pT of W-jet 2
XLabel=$\pT$ (GeV)
#YLabel=$\dfrac{1}{\sigma}\dfrac{\text{d}\sigma}{\text{d}\pT$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/paraproj_asymm_both
Title=Combined parallel asymmetry of pull vector of both $W$-jets
XLabel=$|t_{\parallel}|$
YLabel=$A_{\parallel}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/pull_12_all_highpT
Title=Pull angle 1,2 between $W$-jets, all constituents, jetpT $>$ 60 GeV
RatioPlotYMin=0.95
RatioPlotYMax=1.05
YMin=0.8
XLabel=$\phi(j^W_1, j^W_2) / \pi$
YLabel=$\dfrac{1}{\sigma}\dfrac{\text{d}\sigma}{\text{d}\phi(j^W_1,j^W_2)}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/pull_12_all_veryhighpT
Title=Pull angle 1,2 between $W$-jets, all constituents, jetpT $>$ 80 GeV
RatioPlotYMin=0.95
RatioPlotYMax=1.05
YMin=0.8
XLabel=$\phi(j^W_1, j^W_2) / \pi$
YLabel=$\dfrac{1}{\sigma}\dfrac{\text{d}\sigma}{\text{d}\phi(j^W_1,j^W_2)}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/paraproj_asymm_highpT
Title=Parallel asymmetry of pull vector of highest p_t $W$-jet, jetpT $>$ 60 GeV
RatioPlotYMin=0.8
RatioPlotYMax=1.2
XLabel=$|t_{\parallel}|$
YLabel=$A_{\parallel}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/paraproj_asymm_veryhighpT
Title=Parallel asymmetry of pull vector of highest p_t $W$-jet, jetpT $>$ 80 GeV
RatioPlotYMin=0.8
RatioPlotYMax=1.2
XLabel=$|t_{\parallel}|$
YLabel=$A_{\parallel}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/paraproj_asymm_onebin_jetpTbins_twobins
Title=Parallel asymmetry of pull vector of highest p_t $W$-jet in jet 1 pT bins
XLabel=$\pT$ (GeV)
YLabel=$A_{\parallel}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/paraproj_asymm_j2_onebin_jetpTbins_twobins
Title=Parallel asymmetry of pull vector of $W$-jet 2 in jet 2 pT bins
XLabel=$\pT$ (GeV)
YLabel=$A_{\parallel}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/bparaproj_asymm
Title=Parallel asymmetry of pull vector of highest p_t $B$-jet
XLabel=$|t_{\parallel}|$
YLabel=$A_{\parallel}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/bparaproj_asymm_onebin
Title=Parallel asymmetry of pull vector of highest p_t $B$-jet
XLabel=$|t_{\parallel}|$
YLabel=$A_{\parallel}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/bperpproj_asymm
Title=Perpendicular asymmetry of pull vector of highest p_t $B$-jet
XLabel=$|t_{\bot}|$
YLabel=$A_{\bot}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/bperpproj_asymm_onebin
Title=Perpendicular asymmetry of pull vector of highest p_t $B$-jet
XLabel=$|t_{\bot}|$
YLabel=$A_{\bot}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/bparaproj_asymm_onebin_jetpTbins
Title=Parallel asymmetry of pull vector of highest p_t $B$-jet in jet 1 pT bins
XLabel=$\pT$ (GeV)
YLabel=$A_{\parallel}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/bpull_angle_paraproj_asymm_cor
Title=Correlation between B pull angle (mod $\pi$/2) and parallel asymmetry
XLabel=$\phi(j^B_1, j^B_2) / \pi$
YLabel=$A_{\parallel}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/bJet1_pT
Title=pT of B-jet 1 
XLabel=$\pT$ (GeV)
#YLabel=$\dfrac{1}{\sigma}\dfrac{\text{d}\sigma}{\text{d}\pT$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/bJet2_pT
Title=pT of B-jet 2
XLabel=$\pT$ (GeV)
#YLabel=$\dfrac{1}{\sigma}\dfrac{\text{d}\sigma}{\text{d}\pT$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/paraproj_asymm_onebin_lowpT
Title=Parallel asymmetry of pull vector of highest p_t $W$-jet, jetpT $<$ 60 GeV
XLabel=$|t_{\parallel}|$
YLabel=$A_{\parallel}$
RatioPlotYMin=0.7
RatioPlotYMax=1.3
# END PLOT

# BEGIN PLOT /MC_ColourFlow/paraproj_asymm_onebin_highpT
Title=Parallel asymmetry of pull vector of highest p_t $W$-jet, jetpT $>$ 60 GeV
XLabel=$|t_{\parallel}|$
YLabel=$A_{\parallel}$
RatioPlotYMin=0.7
RatioPlotYMax=1.3
# END PLOT

# BEGIN PLOT /MC_ColourFlow/paraproj_asymm_onebin_veryhighpT
Title=Parallel asymmetry of pull vector of highest p_t $W$-jet, jetpT $>$ 80 GeV
XLabel=$|t_{\parallel}|$
YLabel=$A_{\parallel}$
RatioPlotYMin=0.7
RatioPlotYMax=1.3
# END PLOT

# BEGIN PLOT /MC_ColourFlow/bparaproj_asymm_onebin_M
Title=(weighted) Parallel asymmetry of pull vector of highest p_t $B$-jet
XLabel=$|t_{\parallel}|$
YLabel=$A_{\parallel}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/pull_angle_paraproj_asymm_cor_binone
Title=Parallel asymmetry for $W$-jet 1 for pull vectors with pull angle (mod $\pi$/2) $<$ $\pi$/6
XLabel=$\phi(j^W_1, j^W_2) / \pi$
YLabel=$A_{\parallel}$
RatioPlotYMin=0.8
RatioPlotYMax=1.20
# END PLOT

# BEGIN PLOT /MC_ColourFlow/pullmagnitudej1_all
Title=Magnitude of pull vector of highest p_t $W$-jet
XLabel=$|t|$
YLabel=$\dfrac{1}{\sigma}\dfrac{\text{d}\sigma}{\text{d}|t|}$
RatioPlotYMin=0.95
RatioPlotYMax=1.05
# END PLOT

# BEGIN PLOT /MC_ColourFlow/bpull_both_all
Title=Pull angle (1,2 and 2,1) between $B$-jets, all constituents
RatioPlotYMin=0.95
RatioPlotYMax=1.05
YMin=0.8
XLabel=$\phi(j^b_1, j^b_2) / \pi$ , $\phi(j^b_2, j^b_1) / \pi$
YLabel=$\dfrac{1}{\sigma}\dfrac{\text{d}\sigma}{\text{d}\phi(j^b_1,j^b_2)}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/pullvectorj1_perpproj_all
Title=Magnitude of perpendicular projection of pull vector of highest p_t $W$-jet
XLabel=$|t_{\bot}|$
YLabel=$\dfrac{1}{\sigma}\dfrac{\text{d}\sigma}{\text{d}|t_{\bot}|}$
RatioPlotYMin=0.95
RatioPlotYMax=1.05
# END PLOT

# BEGIN PLOT /MC_ColourFlow/pullvectorj1_paraproj_all
Title=Magnitude of parallel projection of pull vector of highest p_t $W$-jet
XLabel=$|t_{\parallel}|$
YLabel=$\dfrac{1}{\sigma}\dfrac{\text{d}\sigma}{\text{d}|t_{\parallel}|}$
RatioPlotYMin=0.95
RatioPlotYMax=1.05
# END PLOT

# BEGIN PLOT /MC_ColourFlow/pull_both_all_lowpT
Title=Pull angle (1,2 and 2,1) between $W$-jets, all constituents, jetpT $<$ 60 GeV
RatioPlotYMin=0.90
RatioPlotYMax=1.10
YMin=0.8
XLabel=$\phi(j^w_1, j^w_2) / \pi$ , $\phi(j^w_2, j^w_1) / \pi$
YLabel=$\dfrac{1}{\sigma}\dfrac{\text{d}\sigma}{\text{d}\phi(j^w_1,j^w_2)}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/pull_both_all_highpT
Title=Pull angle (1,2 and 2,1) between $W$-jets, all constituents, jetpT $>$ 60 GeV
RatioPlotYMin=0.9
RatioPlotYMax=1.10
YMin=0.8
XLabel=$\phi(j^w_1, j^w_2) / \pi$ , $\phi(j^w_2, j^w_1) / \pi$
YLabel=$\dfrac{1}{\sigma}\dfrac{\text{d}\sigma}{\text{d}\phi(j^w_1,j^w_2)}$
# END PLOT

# BEGIN PLOT /MC_ColourFlow/pull_both_all_veryhighpT
Title=Pull angle (1,2 and 2,1) between $W$-jets, all constituents, jetpT $>$ 80 GeV
RatioPlotYMin=0.90
RatioPlotYMax=1.10
YMin=0.8
XLabel=$\phi(j^w_1, j^w_2) / \pi$ , $\phi(j^w_2, j^w_1) / \pi$
YLabel=$\dfrac{1}{\sigma}\dfrac{\text{d}\sigma}{\text{d}\phi(j^w_1,j^w_2)}$
# END PLOT

## paraproj_asymm_both_onebin
## bpullmagnitudej1_all
## bpull_12_21_cor
## leading_w_b_pull_all