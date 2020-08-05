# BEGIN PLOT /MC_TTBAR/*
LogX=0
LogY=0
NormalizeToIntegral=1
RatioPlotYLabel=MC/default
# END PLOT

# BEGIN PLOT /MC_TTBAR_HADRON/jet_mult
XLabel=Jet multiplicity
# END PLOT

# BEGIN PLOT /MC_TTBAR_HADRON/.*_mass
XLabel=$m$ [GeV]
YLabel=$1/\sigma \, \mathrm{d}\sigma/\mathrm{d}m$ [GeV$^{-1}$]
LogY=0
# END PLOT

# BEGIN PLOT /MC_TTBAR_HADRON/.*_dR
XLabel=$\Delta R$
YLabel=$1/\sigma \, \mathrm{d}\sigma/\mathrm{d}\Delta R$
# END PLOT

# BEGIN PLOT /MC_TTBAR_HADRON/.*_eta
XLabel=$\eta$
YLabel=$1/\sigma \, \mathrm{d}\sigma/\mathrm{d}\eta$
# END PLOT

# BEGIN PLOT /MC_TTBAR_HADRON/.*_deta
XLabel=$\Delta \eta$
YLabel=$1/\sigma \, \mathrm{d}\sigma/\mathrm{d}\Delta \eta$
# END PLOT

# BEGIN PLOT /MC_TTBAR_HADRON/.*_dphi
XLabel=$\Delta \phi$
YLabel=$1/\sigma \, \mathrm{d}\sigma/\mathrm{d}\Delta \phi$
# END PLOT

# BEGIN PLOT /MC_TTBAR_HADRON/jet.?_[1234]_pT
XLabel=$p_T$ [GeV]
YLabel=$1/\sigma \, \mathrm{d}\sigma/\mathrm{d}p_T$ [GeV$^{-1}$]
#LogX=1
# END PLOT

# BEGIN PLOT /MC_TTBAR_HADRON/jet_1_pT
Title=Transverse momentum distribution for jet 1
# END PLOT

# BEGIN PLOT /MC_TTBAR_HADRON/jet_2_pT
Title=Transverse momentum distribution for jet 2
# END PLOT

# BEGIN PLOT /MC_TTBAR_HADRON/jet_3_pT
Title=Transverse momentum distribution for jet 3
# END PLOT

# BEGIN PLOT /MC_TTBAR_HADRON/jet_4_pT
Title=Transverse momentum distribution for jet 4
# END PLOT

# BEGIN PLOT /MC_TTBAR_HADRON/jet_HT
Title=$H_T$ distribution for all jets
XLabel=$H_T$ [GeV]
YLabel=$1/\sigma \, \mathrm{d}\sigma/\mathrm{d}H_T$ [GeV$^{-1}$]
# END PLOT

# BEGIN PLOT /MC_TTBAR_HADRON/jetb_1_pT
Title=Transverse momentum distribution for $b$-jet 1
# END PLOT

# BEGIN PLOT /MC_TTBAR_HADRON/jetb_2_pT
Title=Transverse momentum distribution for $b$-jet 2
# END PLOT

# BEGIN PLOT /MC_TTBAR_HADRON/jetl_1_pT
Title=Transverse momentum distribution for light jet 1
# END PLOT

# BEGIN PLOT /MC_TTBAR_HADRON/jetl_2_pT
Title=Transverse momentum distribution for light jet 2
# END PLOT

# BEGIN PLOT /MC_TTBAR_HADRON/W_mass
Title=Mass distribution for $W$ bosons
XLabel=$m_{jj}$ [GeV]
YLabel=$1/\sigma \, \mathrm{d}\sigma/\mathrm{d}m_{jj}$ [GeV$^{-1}$]
LogY=0
Rebin=1
# END PLOT

# BEGIN PLOT /MC_TTBAR_HADRON/t_mass
Title=Mass distribution for reconstructed top
XLabel=$m_{jjb}$ [GeV]
YLabel=$1/\sigma \, \mathrm{d}\sigma/\mathrm{d}m_{jjb}$ [GeV$^{-1}$]
LogY=0
Rebin=1
# END PLOT

# BEGIN PLOT /MC_TTBAR_HADRON/t_mass_W_cut
Title=Mass distribution for reconstructed top after $m_W$ cut
XLabel=$m_{jjb}$ [GeV]
YLabel=$1/\sigma \, \mathrm{d}\sigma/\mathrm{d}m_{jjb}$ [GeV$^{-1}$]
LogY=0
Rebin=1
# END PLOT

# BEGIN PLOT /MC_TTBAR_HADRON/t_mass_W_scaled
Title=Mass distribution for "fitted" top after $m_W$ cut
XLabel=$m_{Wb}$ [GeV]
YLabel=$1/\sigma \, \mathrm{d}\sigma/\mathrm{d}m_{Wb}$ [GeV$^{-1}$]
Rebin=1
# END PLOT

# BEGIN PLOT /CMS_TOP_13_007/.*
LogY=0
# END PLOT

# BEGIN PLOT /CMS_TOP_13_007/chMult*
Title=Charged multiplicity
XLabel=$N_{ch}$
# END PLOT

# BEGIN PLOT /CMS_TOP_13_007/chSumPt*
Title=Charged transverse momentum flux
XLabel=$\Sigma p_T$ [GeV]
# END PLOT
