from sklearn.model_selection import train_test_split
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import pyarrow.parquet as pq
import math

#random.seed(42)
#np.random.seed(42)
tf.random.set_seed(42)
#tf.keras.utils.set_random_seed(42)


# In[160]:
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--cat", type = str, help="Categories: 1mu1tau, 1ele1tau, 2tau", required=True)

args = parser.parse_args()
cat =  args.cat

PLOT_FOLDER = './HZA_plots_dnn/'+cat+'/'

skim = pq.read_table("./parquet/"+cat+"_looseCuts.parquet")
df1 = skim.to_pandas()

print(df1)
print(df1.columns)

df1["label"] = 0
df1.loc[df1.process.str.contains("HZA_signal_600_100"), ['label']] = 1
df1.loc[df1.process.str.contains("HZA_signal_600_150"), ['label']] = 2
df1.loc[df1.process.str.contains("HZA_signal_600_200"), ['label']] = 3
df1.loc[df1.process.str.contains("HZA_signal_600_250"), ['label']] = 4
df1.loc[df1.process.str.contains("HZA_signal_600_300"), ['label']] = 5
df1.loc[df1.process.str.contains("HZA_signal_600_350"), ['label']] = 6
df1.loc[df1.process.str.contains("HZA_signal_600_400"), ['label']] = 7

#some plots, add more

var = {"Zcand_m" :                   [60., 120., 50 ] ,
       "Zcand_pt":                   [0. , 500., 50 ] ,
       "tau1_Hcand_pt":              [0. , 500., 50 ] ,
       "tau1_Hcand_deepTauVSjet":    [0. , 1.  , 30 ] ,
       "tau2_Hcand_pt":              [0. , 500., 50 ] ,
       "tau2_Hcand_deepTauVSjet":    [0. , 1.  , 30 ] ,
       "Hcand_m":                    [0. , 500., 50] ,
       "ZHcand_m":                   [0. , 700., 50]
       }
import matplotlib.pyplot as plt
for key in var:

    plt.hist(df1.loc[df1.process.str.contains("600_100") == True,:][key], histtype=("step"),range=(var[key][0],var[key][1]),  bins=var[key][2], density = True, label = "Signal 600 100" )

    plt.hist(df1.loc[df1.process.str.contains("HZA") == False,:][key], histtype=("step"), range=(var[key][0],var[key][1]),  bins=var[key][2], density= True, label = "Background" )

    plt.legend(loc='best')
    plt.xlabel(key)
    plt.savefig(PLOT_FOLDER+'h_'+key+'.pdf')
    print("Saved")
    #plt.show()
    plt.clf()


#add more vars 
input_vars = ["met_pt"]


print("INPUT VAR for DNN")
print(input_vars)

print("---------------------------------")
print("events before cuts")
print("Bkg")
b = df1.loc[df1.label ==0,"weight"].sum()
print(b)
print("Signal")
for i in range(1,8):
     print(i)
     s = df1.loc[df1.label ==i,"weight"].sum()
     print(s)
     print("S/B", s/math.sqrt(b))

#SIGNAL REGION

#modify the cuts, put the ones you choose
sel_tightZ = (df1.lep1_Zcand_pfRelIso04_all < 0.15 ) & (df1.lep1_Zcand_tightId) & (df1.lep2_Zcand_pfRelIso04_all < 0.15 ) & (df1.lep2_Zcand_mediumId) 
sel_tau2   = (df1.tau2_Hcand_iddeepTauVSmu >= 1 )  & (df1.tau2_Hcand_iddeepTauVSe >= 1 ) & (df1.tau2_Hcand_iddeepTauVSjet >= 1 )
if cat == "1mu1tau": 
   sel_tau1 = (df1.tau1_Hcand_mediumId ) & (df1.tau1_Hcand_pfRelIso04_all < 0.25)
elif cat == "1ele1tau":
   sel_tau1 = (df1.tau1_Hcand_mvaFall17V2noIso_WP90) #& (df1.tau1_Hcand_pfRelIso03_all < 0.25)
elif cat == "2tau":
   sel_tau1 = (df1.tau1_Hcand_iddeepTauVSmu >= 1 )  & (df1.tau1_Hcand_iddeepTauVSe >= 1 ) &  (df1.tau1_Hcand_iddeepTauVSjet >= 1 )

OS = df1.tau1_Hcand_ch != df1.tau2_Hcand_ch 
selection = sel_tightZ & sel_tau1 & sel_tau2 & OS

df= df1.loc[selection, :]

print("---------------------------------")
print("events after cuts")
print("Bkg")
b = df.loc[df.label ==0,"weight"].sum()
print(b)
print("Signal")
for i in range(1,8):
     print(i)
     s = df.loc[df.label ==i,"weight"].sum()
     print(s)
     print("S/B", s/math.sqrt(b))

print("-----------------------------------------")
print("-----------------------------------------")
print(df.process.unique())
print(df[input_vars])
print(df.label)


# In[48]:


X_train, X_test, y_train, y_test = train_test_split(df[input_vars], df.label,
                                                    test_size=0.30, random_state=42)

#plot correlation matrix between input vars

#choose a dnn

inputs = keras.Input(shape=(len(input_vars),), name="particles")
x = layers.Dense(64, activation="relu", name="dense_1")(inputs)
outputs = layers.Dense(8, activation="softmax", name="predictions")(x)#, activation="softmax", name="predictions")(x)

model = keras.Model(inputs=inputs, outputs=outputs)

model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])

earlystop = tf.keras.callbacks.EarlyStopping(monitor='val_accuracy', min_delta=0.001, patience=10, verbose=1, mode='auto')

Y = tf.keras.utils.to_categorical(y_train)

history = model.fit(
    X_train,
    Y,
    batch_size=128,
    epochs=100,
    # We pass some validation for
    # monitoring validation loss and metrics
    # at the end of each epoch
    callbacks = [earlystop],   
    validation_data=(X_test, y_test),
    validation_split = 0.1
)

model.save(PLOT_FOLDER+"dnn_multiclass")


plt.figure(figsize=(10, 7))
plt.subplot(1,2,1)
plt.plot(history.history['loss'],color='m',label='Training loss')
plt.plot(history.history['val_loss'],color='b',label='Validation loss')
plt.xlabel('Epochs')
plt.ylabel('Loss function')
plt.ylim(ymin = 0.0)
plt.legend()

#plot accurancy

predictions_conf = model.predict(X_test)


print(predictions_conf.shape)
print(y_test.shape)

response, _, _ = np.histogram2d(np.argmax(predictions_conf, axis=1), y_test, bins=[np.linspace(0,8,9),np.linspace(0,8,9)])     

sum_of_rows = response.sum(axis=1)
norm_response = response / sum_of_rows[:, np.newaxis]

fig = plt.figure(figsize=(16, 12))
ax2 = fig.add_subplot(111)
size=8
x_start = 0.0
x_end = 8.0
y_start = 0.0
y_end = 8.0
extent = [x_start, x_end, y_start, y_end]
im=ax2.imshow(norm_response,extent=extent,origin="lower",interpolation='None',cmap="cool")
ax2.set_xlabel("Prediction")
ax2.set_ylabel("True")

class_names = ['bkg', 'HZA_600_100', 'HZA_600_150', 'HZA_600_200', 'HZA_600_250', 'HZA_600_300', 'HZA_600_350', 'HZA_600_400']

ax2.set_xticks(range(8))
ax2.set_yticks(range(8))
ax2.set_xticklabels(class_names,rotation_mode="anchor", ha="left")
ax2.set_yticklabels(class_names,rotation_mode="anchor", ha="left", rotation=90)

# Add the text
jump_x = (x_end - x_start) / (2.0 * size)
jump_y = (y_end - y_start) / (2.0 * size)
x_positions = np.linspace(start=x_start, stop=x_end, num=size, endpoint=False)
y_positions = np.linspace(start=y_start, stop=y_end, num=size, endpoint=False)

for y_index, y in enumerate(y_positions):
    for x_index, x in enumerate(x_positions):
        label = round(norm_response[y_index, x_index],3)
        text_x = x + jump_x
        text_y = y + jump_y
        ax2.text(text_x, text_y, label, color='black', ha='center', va='center')
        
fig.colorbar(im)
fig.savefig(PLOT_FOLDER+'confusion_matrix.pdf')
plt.clf()



plt.hist(1-predictions_conf[:,0],  label='bkg', histtype=("step"), bins =50,density=True);
plt.hist(1-predictions_conf[:,1],  label='signal 1', histtype=("step"), bins =50,density=True);
plt.hist(1-predictions_conf[:,2],  label='signal 2', histtype=("step"), bins =50,density=True);
plt.hist(1-predictions_conf[:,3],  label='signal 3', histtype=("step"), bins =50,density=True);
plt.hist(1-predictions_conf[:,4],  label='signal 4', histtype=("step"), bins =50,density=True);
plt.hist(1-predictions_conf[:,5],  label='signal 5', histtype=("step"), bins =50,density=True);
plt.hist(1-predictions_conf[:,6],  label='signal 6', histtype=("step"), bins =50,density=True);
plt.hist(1-predictions_conf[:,7],  label='signal 7', histtype=("step"), bins =50,density=True);

plt.legend(loc = 'upper center')
plt.xlabel('DNN score')
plt.yscale('log')

plt.yscale('log')
plt.savefig(PLOT_FOLDER+'dnn_score_log.pdf')


#plot ROC curves
