# ENP_Post-Fire_Recovery

Naming convention for scripts: PRODUCT_Description_Author.R

	•	Figure_BaselineModel_SLM – Creates model results figures for the baseline model
	•	Figure_BModelOBSvsPRED_SLM – Observed versus predicted model figures

Data wrangling and analysis codes in “Scripts” folder. 
Naming convention: ORDER # _ Description_ Author.R
1 = data prep
2 = baseline
3 = Recovery
4 = Drivers

	•	1.1_VegLayers_MGM
	•	3.2_RecoveryRates_MGM

## FIGURE FORMATTING:
Text size: 20. 
Theme: BW.    

### NDVI / Thresholds
Greens:  MetBrewer “VanGogh3” : Low  (light green) > High (dark green). 
#e7e5cc,  #c3d6a3, #9bc185, #669e62, #3d7d3d, #1f5b25, #1e3e14, #1a2914

### Fire frequency
Warm colors: MetBrewer “Tam”:  Low (yellow)  > High (purple).    
#fed352, #ffb142, #f08737, #df4f34, #bb292c, #a12d54, #611f5d, #341648

### Model evaluation. 
MetBrewer “Veronese. 
Modeled: #6f948c. 
Observed: #2d6b68. 
Target: #67322e.    
Pre-fire NDVI: #c38f17. 
Post-fire NDVI: #155449. 

### Climate
MetBrewer “Veronese.   
Wet: #112c42 / #738e8e. 
Normal: #6f948c / #677853. 
Dry: #99610a / #a9845a. 

#Sample GGPLOT Code:
```r{
ggplot() +
  …………+
  labs(
    x=“…”,
    y= “…”,
    title = “….”,
    color=“….”, 
    fill= “….”) + 
  #scale_color_met_c(“PALLET NAME “) +
  # scale_fill_met_d(“PALLET NAME”) +
  #scale_fill_manual(values=c(#, #, #) +
  theme_bw()+
  theme(
    text = element_text(size = 20),
    legend.position=“right”)
  )
}```
