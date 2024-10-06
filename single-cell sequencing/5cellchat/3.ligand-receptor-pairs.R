# View cell subgroups and factor order:
levels(cellchat@idents)
# Select the cell subgroups of interest:
vertex.receiver = c(12,13,14,15,16,19,21,22,23,24,25)
par(mfrow = c(1,1))
netVisual_aggregate(cellchat,
                    layout = c('hierarchy'), #"circle", "hierarchy", "chord"
                    signaling = pathways.show,
                    vertex.receiver = vertex.receiver)

pdf(file = "contribution.pdf", width =10 , height = 30)
# Calculate the contribution of ligand-receptor pairs in the target signaling pathway:
netAnalysis_contribution(cellchat, signaling = pathways.show) # Bar chart for ligand-receptor pair contributions
dev.off()
##################################################################################
# Extract ligand-receptor pairs:
pairLR.CXCR4 <- extractEnrichedLR(cellchat,
                                 signaling = pathways.show,
                                 geneLR.return = FALSE)
LR.show <- pairLR.CXCR4[1,] # Taking the top 1 ligand-receptor pair by contribution as an example
pairLR.CXCR4; LR.show
# Hierarchy plot:
pdf(file = "MIF_CD74_CXCR4_netVisuall.pdf", width =12 , height = 8)
netVisual_individual(cellchat,
                     layout = c('hierarchy'),
                     signaling = pathways.show, # Target signaling pathway
                     pairLR.use = LR.show, # Target ligand-receptor pair
                     vertex.receiver = vertex.receiver) # Cell subgroups of interest
dev.off()
##################################################################################
# Extract ligand-receptor pairs:
pairLR.CD44 <- extractEnrichedLR(cellchat,
                                  signaling = pathways.show,
                                  geneLR.return = FALSE)
LR.show <- pairLR.CD44[2,] # Taking the top 1 ligand-receptor pair by contribution as an example
pairLR.CD44; LR.show
# Hierarchy plot:
pdf(file = "MIF_CD74_CD44_netVisuall.pdf", width =12 , height = 8)
netVisual_individual(cellchat,
                     layout = c('hierarchy'),
                     signaling = pathways.show, # Target signaling pathway
                     pairLR.use = LR.show, # Target ligand-receptor pair
                     vertex.receiver = vertex.receiver) # Cell subgroups of interest
dev.off()

##################################################################################
# Extract ligand-receptor pairs:
pairLR.FN1CD44 <- extractEnrichedLR(cellchat,
                                 signaling = pathways.show,
                                 geneLR.return = FALSE)
LR.show <- pairLR.FN1CD44[61,] # Taking the top 1 ligand-receptor pair by contribution as an example
pairLR.FN1CD44; LR.show
# Hierarchy plot:
pdf(file = "FN1_CD44_netVisuall.pdf", width =12 , height = 8)
netVisual_individual(cellchat,
                     layout = c('hierarchy'),
                     signaling = pathways.show, # Target signaling pathway
                     pairLR.use = LR.show, # Target ligand-receptor pair
                     vertex.receiver = vertex.receiver) # Cell subgroups of interest
dev.off()