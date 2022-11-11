import sys
import numpy as np 
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

file = sys.argv[1]
k562_model_predictions = []
k562_model_observations = []
descriptions = []
gene_names = []
with open(file) as f:
	for i,line in enumerate(f.readlines()):
		if line.strip('"').startswith('##'):
			header = np.array(line.strip('"\r\n').split('\t'))
			k562_obs_idx = np.where(header == 'E123')[0][0]
			print(k562_obs_idx)
			continue
		elif not line.strip('"').startswith('#'):
			fields = line.strip('"\r\n').split('\t')
			k562_model_predictions.append(float(fields[4]))
			k562_model_observations.append(float(fields[k562_obs_idx]))
			gene_names.append(fields[1])
			descriptions.append(fields[2])

genesoi = ["PIM1", "SMYD3", "FADS1", "PRKAR2B", "GATA1", "MYC"]
genesoilocs =[]
for geneoi in genesoi:
    genesoilocs.append(np.where(np.array(gene_names) == geneoi)[0][0])
for i in range(len(descriptions)):
    if "hemoglobin subunit" in descriptions[i]:
        genesoi.append(gene_names[i])
        genesoilocs.append(i)

corr = pearsonr(k562_model_observations,k562_model_predictions)

fig,ax = plt.subplots()
ax.scatter(k562_model_observations,k562_model_predictions,color='blue',alpha=.2,s=2)
ax.set_xlabel(r'$log_{10}$(K562 expression level)')
ax.set_ylabel(r'$log_{10}$(K562 predicted expression)')
line_xs = np.linspace(max(min(k562_model_predictions) ,min(k562_model_observations) ), min(max(k562_model_predictions) ,max(k562_model_observations)), 100)
line_ys = 0 + 1 * line_xs
ax.plot(line_xs,line_ys,'r--')
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)
ax.text(1.9, 1.5, "r^2 = " + str(round(corr[0]**2,2)) + "\nn = " + str(len(k562_model_observations)))
for i,g in enumerate(genesoi):
	print(g)
	ax.text(k562_model_observations[genesoilocs[i]], k562_model_predictions[genesoilocs[i]], g, color="maroon",fontweight="demi")
	ax.scatter(k562_model_observations[genesoilocs[i]], k562_model_predictions[genesoilocs[i]],color='maroon',s=2,alpha=.2)
plt.show()
