import pandas as pd
#from statsmodels.formula.api import ols
#from statsmodels.stats.anova import anova_lm
from statsmodels.graphics.factorplots import interaction_plot
import matplotlib.pyplot as plt
from scipy import stats




#datafile="/home/vital/Desktop/Quanti_PeptMatrix/one_peptide_test.csv"

datafile="/home/vital/Desktop/Quanti_PeptMatrix/Quanti_PeptMatrix.csv"

data = pd.read_csv(datafile)



#check interaction plot for one particular peptide (datafile = one_peptide_test.csv)
#fig = interaction_plot(data.drug, data.cell_line, data.intensity)
#plt.show()




#DEGREES OF FREEDOM:
def dof(peptide_df):
	N = len(peptide_df.intensity)
	dof_cell_line = len(peptide_df.cell_line.unique()) - 1
	dof_drug = len(peptide_df.drug.unique()) - 1
	dof_b = dof_cell_line*dof_drug #between
	dof_w = N - (len(peptide_df.cell_line.unique())*len(peptide_df.drug.unique())) #within
	return(dof_cell_line, dof_drug, dof_b, dof_w)



#CALCULATION SUM OF SQUARES:
def ssq(data, grand_mean):
	#SUM of SQUARES FOR EACH FACTOR (Independent variables: cell_line AND drug): 
	ssq_cell_line = sum( [(data[data.cell_line == l].intensity.mean() - grand_mean)**2 for l in data.cell_line] )
	ssq_drug = sum( [(data[data.drug == d].intensity.mean() - grand_mean)**2 for d in data.drug] )
	#TOTAL SUM OF SQUARES:
	ssq_t = sum((data.intensity - grand_mean)**2)
	#SUM OF SQUARES WITHIN (Also called residual)
	lung= peptide_df[peptide_df.cell_line == "Lung"]
	mel= peptide_df[peptide_df.cell_line == "Melanoma"]
	col= peptide_df[peptide_df.cell_line == "Colon"]	
	#For the cell_line subsets, get means for each drug:
	lung_drug_means = [lung[lung.drug == d].intensity.mean() for d in lung.drug]
	mel_drug_means = [mel[mel.drug == d].intensity.mean() for d in mel.drug]
	col_drug_means = [col[col.drug == d].intensity.mean() for d in col.drug]
	#And GET SSQ_w
	ssq_w = sum((mel.intensity - mel_drug_means)**2) + sum((lung.intensity - lung_drug_means)**2) + sum((col.intensity - col_drug_means)**2) 
	#SUM OF SQUARES OF THE INTERACTION (BETWEEN)
	ssq_b = ssq_t-ssq_cell_line-ssq_drug-ssq_w
	return(ssq_cell_line, ssq_drug, ssq_t, ssq_w, ssq_b)




#CALCULATION MEAN SQUARES:
def msq(ssq_cell_line, ssq_drug, ssq_b, ssq_w, dof_cell_line, dof_drug, dof_b, dof_w):
	msq_cell_line = ssq_cell_line / dof_cell_line
	msq_drug = ssq_drug / dof_drug
	msq_cell_line_x_drug = ssq_b / dof_b #INTERACTION
	msq_w = ssq_w/dof_w #WITHIN
	return(msq_cell_line, msq_drug, msq_cell_line_x_drug, msq_w)




#CALCULATION OF F-ratio (F-statistic)
#The F-statistic is simply the mean square for each effect and the interaction divided by the mean square for within (error/residual).
def f_ratio(msq_cell_line, msq_drug, msq_cell_line_x_drug, msq_w):
	f_cell_line = msq_cell_line/msq_w
	f_drug = msq_drug/msq_w
	f_cell_line_x_drug = msq_cell_line_x_drug/msq_w
	return(f_cell_line, f_drug, f_cell_line_x_drug)









for index, row in data.iterrows():
	#Get pandas df for each peptide with columns intensity, cell_line and drug:
	print index
	peptide_df = pd.DataFrame(columns=["intensity", "cell_line", "drug"])
	peptide_df.cell_line = ["Melanoma"] * 24 + ["Lung"] * 24 + ["Colon"] * 24
	peptide_df.drug = (["CT0"]*3 + ["5-FU"]*3 + ["DOXO"]*3 + ["MTX"]*3 + ["CT72"]*3 + ["PCTL"]*3 + ["SEN"]*3 + ["TDX"]*3)*3
	peptide_df.intensity = row.values[1:-1]
	#get Degrees of Freedom
	dof_cell_line, dof_drug, dof_b, dof_w = dof(peptide_df)
	#get grand mean Using Pandas DataFrame method mean on the dependent variable only
	grand_mean = peptide_df.intensity.mean() #which is simply the mean of all intensity values
	#Calculate sum of squares
	ssq_cell_line, ssq_drug, ssq_t, ssq_w, ssq_b = ssq(peptide_df, grand_mean)
	#Calculate mean squares
	msq_cell_line, msq_drug, msq_cell_line_x_drug, msq_w = msq(ssq_cell_line, ssq_drug, ssq_b, ssq_w, dof_cell_line, dof_drug, dof_b, dof_w)
	#Calculate F-ratio (F-satistic)
	f_cell_line, f_drug, f_cell_line_x_drug = f_ratio(msq_cell_line, msq_drug, msq_cell_line_x_drug, msq_w)
	#OBTAIN p-values
	#scipy.stats method f.sf checks if the obtained F-ratios are above the critical value. 
	#Use F-value for each effect and interaction as well as the degrees of freedom for them, and the degree of freedom within.
	#Null Hypothesis 1:  H0: cell-lines have the same response
	p_cell_line = stats.f.sf(f_cell_line, dof_cell_line, dof_w)
	#Null Hypothesis 2:  H0: drugs have the same response
	p_drug = stats.f.sf(f_drug, dof_drug, dof_w)
	#Null Hypothesis 3:  H0: cell-lines do not interact with drugs in the response
	p_cell_line_x_drug = stats.f.sf(f_cell_line_x_drug, dof_b, dof_w)
	#ADD rows to the original pandas df with the three pvalues:
	data.loc[index, 'p_cell_line'] = p_cell_line
	data.loc[index, 'p_drug'] = p_drug
	data.loc[index, 'p_cell_line_x_drug'] = p_cell_line_x_drug
	

data.to_csv("Quanti_PeptMatrix_p.csv")


