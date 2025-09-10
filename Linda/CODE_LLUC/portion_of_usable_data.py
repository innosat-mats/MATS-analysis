#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from mats_utils.geolocation.coordinates import heights, fast_heights
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT
import time
from shapely.geometry import Polygon
from tqdm import tqdm

import pandas as pd
import plotly.express as px
#%%
# This generates the "filter" dictionary that can be used to extract images of the specified channel only
def get_filter(channel):
    
    """
    Function that select the channel that is going to be used
    Args: - channel: channel name
    Returns: dict
    """


    filters = {"IR1": 1, "IR2": 4, "IR3": 3, "IR4": 2, "UV1": 5, "UV2": 6}
    try:
        filt = filters[channel]
    except Exception:
        raise ValueError(f"Invalid channel: {channel}!")

    return {'CCDSEL': [filt, filt]}

# # This extracts the specified variables for one image (image number idx) from pandas data frame 
# # and packs them into a simple python dictionary (no Pandas any more!)
# def get_image(data, idx, var):
#     """
#     Function that returns data of a single image with the specified variables
#     Args: 
#     - data: df with necessary data
#     - idx: number of image
#     - list of variables given
#     Returns: list values of variables of a single image
#     """
    
#     res = {"num_image": idx}
#     for v in var:
#         if len(data[v].shape) > 1:
#             res[v] = data[v][idx, ...]
#         else:
#             res[v] = data[v][idx]
#     return res

#Function that calculates the proportion of good heights for a given position


def interpolar_altura(matriz, punto):

    """
    Function that interpolates heights values between two points (two pixels)
    Args: 
    - matriz: array of heights
    - punto: point going to be interpolated (it can be a decimal number)
    Returns:
    - altura_interpolada: value of interpolated height
    """

    # Get the integer coordinates of the neighbouring points
    x0, y0 = int(punto[0]), int(punto[1])
    x1, y1 = x0 + 1, y0 + 1
    
    # Obtain the height values of neighbouring points
    z00, z01 = matriz[y0, x0], matriz[y0, x1]
    z10, z11 = matriz[y1, x0], matriz[y1, x1]
    
    # Linearly interpolate in the x-direction
    z0 = z00 + (z01 - z00) * (punto[0] - x0)
    z1 = z10 + (z11 - z10) * (punto[0] - x0)
    
    # Linearly interpolate in the y-direction
    altura_interpolada = z0 + (z1 - z0) * (punto[1] - y0)
    
    return altura_interpolada   


def encontrar_borde_con_gradiente(matriz, punto, gradiente):

    """
    Function that finds the edge or the limit value of 50km from a given point and gradient
    Args:
    - matriz: array of heights
    - punto: point where it starts from
    - gradiente: gradient in "punto"
    Returns:
    - intersection point
    - track: track followed to find the edge with pixel coordinates
    - heights_track: track of heights followed to find the intersection
    """

    y, x = punto
    
    # To normalise the vector of the gradient
    gradiente = (gradiente[0]/2/2/2, gradiente[1]/40/40)

    norma = np.sqrt((gradiente[0])**2+(gradiente[1])**2)

    gradiente_normalizado = gradiente / norma #np.linalg.norm(gradiente)

    
    

    track = [(y,x)]

    heights_track = [heights[punto]]
    
    # Loop until find the intersection
    contador_pasos = 0
    while 0 <= y < matriz.shape[0] and 0 <= x < matriz.shape[1]:

        y += gradiente_normalizado[0]
        x += gradiente_normalizado[1]
        contador_pasos += 1

        point_track = (y,x)

        track.append(point_track)

        if 0 <= y < (matriz.shape[0]-1) and 0 <= x < (matriz.shape[1]-1):
            interpolated_height = interpolar_altura(matriz, (point_track[1], point_track[0]))

            heights_track.append(interpolated_height) 
        
        # Here verifies if the point is out of bounds
        if not (0 <= y < (matriz.shape[0]-1) and 0 <= x <(matriz.shape[1]-1) and matriz[int(round(y)), int(round(x))]>50000):
            # Round the point to be in the array
            y = np.clip(y, 0, matriz.shape[0] - 1)
            x = np.clip(x, 0, matriz.shape[1] - 1)

            track.pop()
            heights_track.pop()
    
            return (int(round(y)), int(round(x))), track , heights_track
        
    
    return None



def calcular_area_poligono(vertices):
    """
    Function that calculates the area of useful data in the image
    Args:
    - vertices: list of starting and intersection points that form a polygon
                (they must be in order! However, the script is done to be like that)
    Returns:
    - area: area of polygon (in pixels**2)
    """
    
    if len(vertices)>2:
        poligono = Polygon(vertices)
        
        # Calcular el área del polígono
        area = poligono.area
    
    else:
        area = 0
    
    return area



def good_data(heights, plot=False):

    """
    Function that calculates the proportion of useful data in a single image
    and the vertical profiles of useful values
    Args:
    - heights: array of heights for a single image
    Returns:
    - area: pixels**2 of useful data
    - proportion of useful values for an image
    - vertical_profiles: array with the vertical profiles of useful data

    IMPORTANT: at the end of the function there is a silent section to plot the image with the useful data and the vertical profiles
    """

    altura_isohipsa = 100000  # en metros

    puntos_isohipsa = []

    # Traverse the height matrix to find the points that match the desired height.
    for i in range(len(heights)):
        for j in range(len(heights[0])):
            if abs(heights[i][j] - altura_isohipsa) < 135:  # Consider as coincidence if the difference is less than 135 metres.
                puntos_isohipsa.append((i, j)) 


    grad_y, grad_x = np.gradient(heights)
    vertical_profiles = []
    intersection = []
    max_height_valid = []

    #here finds all intersection points

    for p1 in puntos_isohipsa:

        p1 = (p1[1], p1[0])
    
        punto, track, heights_track = encontrar_borde_con_gradiente(heights, (p1[1],p1[0]), (-grad_y[p1[1],p1[0]], -grad_x[p1[1],p1[0]]))

        if (len(heights_track)) != 0:
            if (heights_track[-1]<=80000):

                vertical_profiles.append(heights_track)
                intersection.append(punto)
                max_height_valid.append((p1[1], p1[0]))

        if ((len(heights_track)) == 0 or (heights_track[-1]>80000)) :
            vertical_profiles.append([])



    intersection_reverse = []

    for i in range (1,len(intersection)+1):
        intersection_reverse.append(intersection[-i])


    #polygon of useful data
    points_to_plot = max_height_valid+intersection_reverse


    area = calcular_area_poligono(points_to_plot)

    num_rows, num_columns = heights.shape

    proportion = area / (num_rows*num_columns)

    #TO PLOT THE USEFUL DATA IN THE IMAGE AND THE VERTICAL PROFILES
    
    if len(vertical_profiles)!=0:

        longitud_maxima = max(len(lista) for lista in vertical_profiles)

        # Rellena las listas más cortas con NaN hasta que todas tengan la misma longitud
        lista_de_listas_con_nan = [lista + [np.nan] * (longitud_maxima - len(lista)) for lista in vertical_profiles]

        # Convierte la lista de listas en una matriz bidimensional de NumPy
        new_vertical_profiles = np.array(lista_de_listas_con_nan)

    else:
        new_vertical_profiles = vertical_profiles

    if plot:
        plot_contour_with_points(heights, points_to_plot)

        for i, columna in enumerate(new_vertical_profiles, start=1):
            plt.scatter([i] * len(columna), columna, label=f'Columna {i}')

        # Configurar etiquetas y título
        plt.xlim(0.5, len(new_vertical_profiles)+0.5)
        plt.xlabel('Column')
        plt.ylabel('Height (m)')
        plt.title('Vertical profiles')
        #plt.xticks(range(1, len(vertical_profiles) + 1))
        #plt.legend()
        plt.grid(True)

        # Mostrar el gráfico
        plt.show()
    

    return area, proportion, vertical_profiles


    

def list_of_proportions(df):

    """
    Function that creates a list of proportions of useful data with a df
    with the desired variables for a determined period of time
    Args:
    - df: df with data of a determined period of time
    Returns:
    - list of proportions of useful data in an image for all iterations
    """



    proportions = []
    total_iterations = len(df) - 1
    with tqdm(total=total_iterations, desc="Processing") as pbar:
        for i in range(total_iterations):
            image = get_image(df, i, image_vars) # Get all the needed variables
            heights = fast_heights(image)
            area, proportion, vertical_profiles = good_data(heights)
            proportions.append(proportion)
            pbar.update(1)
            remaining_iterations = total_iterations - i - 1
            pbar.set_postfix({"Remaining iterations": remaining_iterations})
            time.sleep(1)
    return proportions



def plot_contour_with_points(alturas, puntos):
    
    """
    Function that plots the image with the useful data
    Args:
    - alturas: array of heights (image)
    - puntos: points to be plotted
    Returns: plot
    """


    num_filas, num_columnas = alturas.shape

    figura_tamano = (num_columnas*0.030*40, num_filas*0.030*2)  # It's been modified due to the size of pixels
    
    plt.figure(figsize=figura_tamano)
    plt.contourf(alturas, cmap='viridis') 
    

    for punto in puntos:
        plt.plot(punto[1], punto[0], 'ro', markersize = 30)  # El índice de la fila es el segundo elemento y el índice de la columna es el primer elemento


    plt.xlabel('Column')
    plt.ylabel('Row')
    plt.title('Useful data')
    plt.grid(True) 

    plt.colorbar()
    
    plt.show()


def guardar_datos(datos, nombre_archivo):
    """
    Function that saves the data as nombre_archivo.npz
    Args:
    - datos: data to be saved
    - nombre_archivo: file name
    """
    np.savez(nombre_archivo, datos=datos)


#%%



"""Specify date, version and channel to be used"""

CHANNEL='IR1'
VERSION='1.0'
START_TIME = [2025, 1, 16, 22, 0, 0]
STOP_TIME = [2025, 1, 17, 2, 0, 0]

for idate in range(23, 29):
    START_TIME1 = [2025, 1, idate, 10, 0, 0]
    STOP_TIME1 = [2025, 1, idate, 12, 0, 0]
    START_TIME2 = [2025, 1, idate+1, 4, 0, 0]
    STOP_TIME2 = [2025, 1, idate+1, 6, 0, 0]
    if idate == 24:
        START_TIME1 = [2025, 1, idate, 10, 0, 0]
        STOP_TIME1 = [2025, 1, idate, 11, 0, 0]
        START_TIME2 = [2025, 1, idate, 11, 0, 0]
        STOP_TIME2 = [2025, 1, idate, 12, 0, 0]
    if idate == 25:
        START_TIME1 = [2025, 1, idate, 13, 0, 0]
        STOP_TIME1 = [2025, 1, idate, 15, 0, 0]
        START_TIME2 = [2025, 1, idate+1, 4, 0, 0]
        STOP_TIME2 = [2025, 1, idate+1, 6, 0, 0]

    dflong1 = read_MATS_data(DT.datetime(*START_TIME1), DT.datetime(*STOP_TIME1),
                        get_filter(CHANNEL), level="1b", version=VERSION)
    dflong2 = read_MATS_data(DT.datetime(*START_TIME2), DT.datetime(*STOP_TIME2),
                        get_filter(CHANNEL), level="1b", version=VERSION)
    #concatenate the two dataframes
    dflong = pd.concat([dflong1, dflong2])
    

    # set df to every 10th row
    df = dflong.iloc[::10]

    # image = get_image(df, 0, image_vars) # Getting all the required variables for the first image (idx=0)
    # heights = fast_heights(image)
    # area, proportion, vertical_profiles=good_data(heights, plot=True)

    #image = get_image(df, 0, image_vars) # Getting all the required variables for the first image (idx=0)
    #
    heights = fast_heights(df.iloc[0]) # This is to avoid a bug which seems related to that heights is a function and a numpy

    # # heights for the first image
    df['imheights'] = df.apply(lambda x: fast_heights(x), axis=1) # Getting heights for all images

    df['area'], df['proportion'], df['vertical_profiles'] = zip(*df['imheights'].apply(good_data)) # Getting area, proportion and vertical profiles for all images

    df.to_pickle('/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/df_portion_jan_'+str(idate)+'.pkl')




    # Define the grid size (in degrees)
    grid_size = 1.0

    # Create a new column for the grid cell each point falls into
    df['lon_grid'] = (df['TPlon'] // grid_size) * grid_size
    df['lat_grid'] = (df['TPlat'] // grid_size) * grid_size

    # Calculate the average proportion for each grid cell
    avg_proportion_grid = df.groupby(['lon_grid', 'lat_grid'])['proportion'].mean().reset_index()

    # Create an interactive map with Plotly
    fig = px.scatter_geo(avg_proportion_grid, lat='lat_grid', lon='lon_grid', 
                        color='proportion', color_continuous_scale='viridis',
                        labels={'proportion': 'Proportion at good height'},
                        title=str(df['TMHeaderTime'].iloc[0]) + ' to ' + str(df['TMHeaderTime'].iloc[-1]))
    fig.show()
    # save the figure as png
    fig.write_image('/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/portion_jan_'+str(idate)+'.png')

# %%
