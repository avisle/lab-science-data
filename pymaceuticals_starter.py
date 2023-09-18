#!/usr/bin/env python
# coding: utf-8

# # Pymaceuticals Inc.
# ---
# 
# ### Analysis
# 
# - Add your analysis here.
#  

# In[15]:


# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st

# Study data files
mouse_metadata_path = "data/Mouse_metadata.csv"
study_results_path = "data/Study_results.csv"

# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)

# Combine the data into a single DataFrame
science_data = pd.merge(study_results, mouse_metadata, how="left", on="Mouse ID")

# Display the data table for preview
science_data


# In[16]:


# Checking the number of mice.
len(science_data["Mouse ID"].unique())


# In[17]:


# Our data should be uniquely identified by Mouse ID and Timepoint
# Get the duplicate mice by ID number that shows up for Mouse ID and Timepoint. 
duplicate_mouse_id = science_data[science_data.duplicated(subset = ["Mouse ID", "Timepoint"])]["Mouse ID"].unique()
duplicate_mouse_id


# In[18]:


# Optional: Get all the data for the duplicate mouse ID. 
duplicated_mouse_data = science_data[science_data["Mouse ID"] == "g989"]
duplicated_mouse_data


# In[19]:


# Create a clean DataFrame by dropping the duplicate mouse by its ID.
cleaned_science_data = science_data[science_data["Mouse ID"].isin(duplicate_mouse_id) == False]
cleaned_science_data


# In[20]:


# Checking the number of mice in the clean DataFrame.
len(cleaned_science_data["Mouse ID"].unique())


# ## Summary Statistics

# In[21]:


# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen

# Use groupby and summary statistical methods to calculate the following properties of each drug regimen: 
# mean, median, variance, standard deviation, and SEM of the tumor volume. 
# Assemble the resulting series into a single summary DataFrame.
means = cleaned_science_data.groupby("Drug Regimen").mean()["Tumor Volume (mm3)"]
median = cleaned_science_data.groupby("Drug Regimen").median()["Tumor Volume (mm3)"]
variance = cleaned_science_data.groupby("Drug Regimen").var()["Tumor Volume (mm3)"]
standard_deviation = cleaned_science_data.groupby("Drug Regimen").std()["Tumor Volume (mm3)"]
SEM = cleaned_science_data.groupby("Drug Regimen").sem()["Tumor Volume (mm3)"]

summary = pd.DataFrame({
    "Mean Tumor Volume": means,
    "Median Tumor Volume": median,
    "Tumor Volume Variance": variance,
    "Tumor Volume Std.Dev.": standard_deviation,
    "Tumor Volume Std. Err.": SEM
})

summary


# In[22]:


# A more advanced method to generate a summary statistics table of mean, median, variance, standard deviation,
# and SEM of the tumor volume for each regimen (only one method is required in the solution)

# Using the aggregation method, produce the same summary statistics in a single line
summary_table = cleaned_science_data.groupby("Drug Regimen").agg({
    "Tumor Volume (mm3)": ["mean", "median", "var", "std", "sem"]
})
summary_table


# ## Bar and Pie Charts

# In[23]:


# Generate a bar plot showing the total number of rows (Mouse ID/Timepoints) for each drug regimen using Pandas.
count = cleaned_science_data["Drug Regimen"].value_counts()
count.plot(kind="bar")
plt.xlabel("Drug Regimen")
plt.ylabel("# of Observed Mouse Timepoints")
plt.xticks(rotation = 90)


# In[24]:


# Generate a bar plot showing the total number of rows (Mouse ID/Timepoints) for each drug regimen using pyplot.
import matplotlib.pyplot as plt 

plt.bar(count.index, count.values)
plt.xlabel("Drug Regimen")
plt.ylabel("# of Observed Mouse Timepoints")
plt.title("Total Number of Rows (Mouse ID/Timepoints) for Each Drug Regimen")
plt.xticks(rotation=90)
plt.tight_layout()

plt.show()


# In[25]:


# Generate a pie plot showing the distribution of female versus male mice using Pandas

gender_distribution = cleaned_science_data["Sex"].value_counts()

gender_distribution.plot(kind="pie", autopct='%1.1f%%')
plt.title("Distribution of Female vs. Male Mice")
 


plt.show()


# In[26]:


# Generate a pie plot showing the distribution of female versus male mice using pyplot
import matplotlib.pyplot as plt
gender_distribution = cleaned_science_data["Sex"].value_counts()

plt.figure(figsize=(6, 6))
plt.pie(gender_distribution, labels=gender_distribution.index, autopct='%1.1f%%')
plt.title("Distribution of Female vs. Male Mice")

plt.show()


# ## Quartiles, Outliers and Boxplots

# In[27]:


# Calculate the final tumor volume of each mouse across four of the treatment regimens:  
# Capomulin, Ramicane, Infubinol, and Ceftamin

# Start by getting the last (greatest) timepoint for each mouse
max_tumor = cleaned_science_data.groupby(["Mouse ID"])["Timepoint"].max()
max_tumor = max_tumor.reset_index()


# Merge this group df with the original DataFrame to get the tumor volume at the last timepoint
merged_data = max_tumor.merge(cleaned_science_data, on=["Mouse ID", "Timepoint"], how="left")


# In[35]:


# Put treatments into a list for for loop (and later for plot labels)
treatments = ["Capomulin", "Ramicane", "Infubinol", "Ceftamin"]

# Create empty list to fill with tumor vol data (for plotting)
tumor_vol_data = []

# Calculate the IQR and quantitatively determine if there are any potential outliers. 
for drug in treatments:
    
    
    # Locate the rows which contain mice on each drug and get the tumor volumes
    final_tumor_vol = merged_data.loc[merged_data["Drug Regimen"] == drug, "Tumor Volume (mm3)"]
    
    # add subset 
    tumor_vol_data.append(final_tumor_vol)
    
    # Determine outliers using upper and lower bounds
    quartiles = final_tumor_vol.quantile([0.25,0.5,0.75])
    lower = quartiles[0.25]
    upper = quartiles[0.75]
    iqr = upper - lower
    lower_bound = lower - (1.5*iqr)
    upper_bound = upper + (1.5*iqr)
outliers = final_tumor_vol.loc[(final_tumor_vol < lower) | (final_tumor_vol > upper)]

print(f"Drug: {drug}")
print(f"Lower Quartile: {lower}")
print(f"Upper Quartile: {upper}")
print(f"IQR: {iqr}")
print(f"Lower Bound: {lower_bound}")
print(f"Upper Bound: {upper_bound}")
print(f"Potential Outliers: {outliers}\n")


# In[36]:


# Generate a box plot that shows the distrubution of the tumor volume for each treatment group.
tumor_vol_data = [merged_data.loc[merged_data["Drug Regimen"] == treatment, "Tumor Volume (mm3)"] for treatment in treatments]
plt.figure(figsize=(10, 6))
plt.boxplot(tumor_vol_data, labels=treatments, sym='r')  # 'sym' adds red outliers
plt.title("Distribution of Tumor Volume by Treatment Group")
plt.xlabel("Treatment Group")
plt.ylabel("Final Tumor Volume (mm3)")
plt.xticks(rotation=45)
plt.grid(axis='y')
plt.show()


# ## Line and Scatter Plots

# In[49]:


# Generate a line plot of tumor volume vs. time point for a single mouse treated with Capomulin
capomulin_data = cleaned_science_data[cleaned_science_data["Drug Regimen"] == "Capomulin"]
mouse_data = capomulin_data[capomulin_data["Mouse ID"] == "l509"]
plt.plot(mouse_data["Timepoint"], mouse_data["Tumor Volume (mm3)"])
plt.xlabel("Timepoint (days)")
plt.ylabel("Tumor Volume (mm3)")
plt.title("Capomulin treatment of mouse l509")


# In[50]:


# Generate a scatter plot of mouse weight vs. the average observed tumor volume for the entire Capomulin regimen
capomulin_data = cleaned_science_data[cleaned_science_data["Drug Regimen"] == "Capomulin"]
capomulin_average = capomulin_data.groupby(["Mouse ID"]).mean()
plt.scatter(capomulin_average["Weight (g)"], capomulin_average["Tumor Volume (mm3)"])
plt.xlabel("Weight (g)")
plt.ylabel("Average Tumor Volume (mm3)")

plt.show()


# ## Correlation and Regression

# In[53]:


# Calculate the correlation coefficient and a linear regression model 
# for mouse weight and average observed tumor volume for the entire Capomulin regimen
correlation = st.pearsonr(capomulin_average["Weight (g)"], capomulin_average["Tumor Volume (mm3)"])
print(f"The correlation between weight and average tumor volume is {round(correlation[0],2)}")

model = st.linregress(capomulin_average["Weight (g)"], capomulin_average["Tumor Volume (mm3)"])
slope = model[0]
b = model [1]
y_values = capomulin_average["Weight (g)"] * slope + b
plt.scatter(capomulin_average["Weight (g)"], capomulin_average["Tumor Volume (mm3)"])
plt.plot(capomulin_average["Weight (g)"], y_values, color = "red")
plt.xlabel("Weight (g)")
plt.ylabel("Average Tumor Volume (mm3)")

plt.show()


# In[ ]:




