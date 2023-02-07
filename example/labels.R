population_labels = c("TBAC" = "A\u03b4-LTMR + A\u03b2-RA-LTMR",
                      "CRTH" = "C-LTMR", "MRTD" = "NP",
                      "CGRT" = "PEP", "TDNV" = "Nociceptors", 
                      "B10.D2" ="b10d2", 
                      "BALB.c" = "balb", 
                      "DRG" = "                                DRG",
                      "iPSC" = "                               iPSC", 
                      "iPSDN" = "iPSDN", "skin" = "skin", 
                      "Aβ.Field.LTMR" = "A\u03b2-Field-LTMR",
                      "Aβ.RA.LTMR" = "A\u03b2-RA-LTMR",
                      "Aβ.SA1.LTMR" = "A\u03b2-SA1-LTMR",
                      "Aδ.LTMR" = "A\u03b4-LTMR",
                      "Nonpeptidergic.Nociceptor"= "NP",
                      "Peptidergic.Nociceptor"= "PEP",
                      "Proprioceptor"         = "PROP",
                      "c"  = "DRG")
sexlabels = c("F" = "Female", "M" = "Male", "mixed" = "Mixed","Mixed" = "Mixed")
subpopulations = c(
  "Nociceptors 3D" = "TDNV_3D.csv", "Nociceptors 4W" = "TDNV_4W.csv",
  "PEP 3D" = "CGRT_3D.csv", "PEP 4W" = "CGRT_4W.csv", 
  "NP 3D" = "MRTD_3D.csv","NP 4W" = "MRTD_4W.csv",
  "C-LTMR 3D" = "CRTH_3D.csv","C-LTMR 4W" = "CRTH_4W.csv",
  "Ad- AB-RA LTMRs 3D" = "TBAC_3D.csv","Ad- AB-RA LTMRs 4W" = "TBAC_4W.csv"
)

subpopulation_labels = c(
  "TDNV_3D" = "Nociceptors 3D", "TDNV_4W" = "Nociceptors 4W",
  "CGRT_3D" = "PEP 3D","CGRT_4W" = "PEP 4W", 
  "MRTD_3D" = "NP 3D", "MRTD_4W" = "NP 4W",
  "CRTH_3D" = "C-LTMR 3D","CRTH_4W" = "C-LTMR 4W",
  "TBAC_3D" = "Ad- AB-RA LTMRs 3D", "TBAC_4W" = "Ad- AB-RA LTMRs 4W", 
  "B10.D2" ="b10d2", 
  "BALB.c" = "balb",
  "DRG" = "                         DRG",
  "old" = "iPSCDN_old",
  "young" = "         iPSCDN_young",
  "Diabetes" = "diabetes",
  "Diabetes_female" = "Diabetes_female",
  "Diabetes_male" = "Diabetes_male",
  "skin" = "skin"
)

# themes 
theme_line = theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=12),
        axis.text.x = element_text(size=10, angle = 45, hjust= 1),
        axis.text.y = element_text(size=10),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(), legend.justification = c(0,0.3)) 

th = theme_bw() + theme(panel.grid = element_blank(),
                        axis.title = element_blank(),
                        axis.text.x = element_text(size=8, angle = 45, hjust= 1),
                        axis.text.y = element_text(size=10),
                        axis.ticks.x = element_blank(),
                        axis.ticks.y = element_blank(), legend.justification = c(0,0.3),
                        strip.text.x = element_text(size = 9, margin = margin(6, 2, 6, 2)), 
                        panel.spacing = unit(0.25, "mm"))


thc = theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 14, vjust=0),
                         plot.subtitle = element_text(size = 12, vjust=0, face="italic"),
                         panel.grid = element_blank(), strip.text = element_text(size = 12),
                         axis.title = element_blank(),
                         axis.text.x = element_text(size=10, angle = 45, hjust= 1),
                         axis.text.y = element_text(size=10),
                         axis.ticks.x = element_blank(),
                         axis.ticks.y = element_blank(), legend.justification = c(0,0.3))