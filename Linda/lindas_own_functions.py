
"""
Created on Fri Dec 2 9:38:36 2022

@author: lindamegner

File for differet plotting functions to be used for MATS commisioning data analysis
"""


import matplotlib.pyplot as plt
import numpy as np

def collapsandplot(imagecube,collapsdim, ax, signallabel='', title='', altvec=None):
    """    
    
    Parameters
    ----------
    imagecub : 3d numpy array made by several CCDimages, created by  create_imagecube(CCDitems, calibrated=False)
    calibrated : BOOLEAN 

    """

    img_hmean=imagecube.mean(collapsdim)
    img_mean=img_hmean.mean(0)
    if altvec is not None:
        myindex1=altvec
    else:
        myindex1=np.arange(0, img_mean.shape[0])
    

    if collapsdim==2:
        myindex=myindex1#*CCDitem['NRBIN']#/7.6+54
        ax.plot(img_mean,myindex, label='mean')

        img_std=img_hmean.std(0)
        ax.plot(img_mean+img_std,myindex, '--', label='mean+1std')
        ax.plot(img_mean-img_std,myindex, '--', label='mean+1std')


        ax.plot(img_hmean.max(0),myindex, '.', label='max')
        ax.plot(img_hmean.min(0),myindex, '.', label='min')

        ax.set_xlabel(signallabel)
        ax.set_ylabel('bin nr')
        ax.set_title(title)
        ax.legend()
        plt.tight_layout()
    elif collapsdim==1:
        myindex=myindex1#*CCDitem['NCBINCCDColumns']#/7.6-125
        ax.plot(myindex, img_mean, label='mean')

        img_std=img_hmean.std(0)
        ax.plot(myindex,img_mean+img_std, '--', label='mean+1std' )
        ax.plot(myindex,img_mean-img_std, '--', label='mean-1std')


        ax.plot(myindex,img_hmean.max(0), '.', label='max')
        ax.plot(myindex,img_hmean.min(0), '.', label='min')

        ax.set_ylabel(signallabel)
        ax.set_xlabel('bin nr')
        ax.set_title(title)
        
        ax.legend()
        plt.tight_layout()

    else:
        raise Warning('collapsdim must be 1 or 2')

    return

def create_imagecube(CCDitems, image_specification='IMAGE'):
    """    
    Parameters
    ----------
    CCDitems : LIST of CCDitems
    image_specification : STR IMAGE or image_calibrated supported 

    Returns
    -------
    3d ndarray with all images and time as the last dimension

    """
    if image_specification!='IMAGE' and image_specification!='image_calibrated':
        Warning('image_specification must be "IMAGE" or "image-calibrated"')
    imagelist=[]
    #for CCDitem in CCDitems:
    for index, CCDitem in CCDitems.iterrows():    
        image=CCDitem[image_specification]
        imagelist.append(image)

    imagecube=np.array(imagelist)

    return imagecube

##########
    
def select_CCDitems(CCDitems, key, value): 
    """
    Function to seache through image items to find specific settings

    Parameters
    ----------
    CCDitems : list
        DESCRIPTION.
    key : string: name og CCDsetting or quanitity, i.e. CCDitem dictonary intex, e.g. 'channel'
        DESCRIPTION.
    value : Value of CCDsetting or quanitity, e.g. 'IR4'
        DESCRIPTION.

    Returns
    -------
    CCDitems_select: List of CCDitems

    """

    CCDitems_select= list(filter(lambda x: ( x[key]==value),CCDitems))
    return CCDitems_select

def select_CCDitems_using_list(CCDitems, key, valuelist):

    
    for value in valuelist:
        CCDitems=select_CCDitems(CCDitems, key, value)
    
    return CCDitems

def rename_CCDitem_entries(df):
    df['NCBIN CCDColumns'] = df['NCBINCCDColumns']
    df['GAIN Truncation'] = df['GAINTruncation']
    df['NCBIN FPGAColumns'] = df['NCBINFPGAColumns']
    mylist=[]
    df['BC'] = [mylist for i in df.index]
    df['GAIN Mode'] = df['GAINMode'] 
    df['EXP Date'] = df['EXPDate'] 
    return df

