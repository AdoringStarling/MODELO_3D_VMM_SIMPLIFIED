import dash
import dash_bootstrap_components as dbc
from dash import html
from dash import dcc
import pandas as pd
import plotly.graph_objects as go
import numpy as np
import os
from plotly.subplots import make_subplots
import io

from geoseismo import *
#Area del estudio
los=-73
loi=-74.4
lai=6.5
las=9

#Se cargan los datos de elevacion estos fueron descargados en https://portal.opentopography.org/datasets
# Global Bathymetry and Topography at 15 Arc Sec: SRTM15+ V2.1  
df_topo   =pd.read_csv('datasets/topo_src_15.xyz',delim_whitespace=True,header=None,decimal='.')
df_topo   =df_topo[(df_topo[1]>lai)&(df_topo[1]<las)&(df_topo[0]>loi)&(df_topo[0]<los)] #Filtros previos
mesh_topo = (df_topo.pivot(index=1, columns=0,values=2))
z_topo,x_topo,y_topo=mesh_topo.values,mesh_topo.columns,mesh_topo.index

topog=go.Surface(z=z_topo,showscale=False, x=x_topo, y=y_topo,colorscale=['black','black'],lighting=dict(ambient=0.3,diffuse=0.5),
                    showlegend=False,opacity=1,name='Topografía',hoverinfo='none')
#Base de datos de sismos convertidos a csv desde http://bdrsnc.sgc.gov.co/paginas1/catalogo/Consulta_Valle_Medio/valle_medio.php
#df_sismos=pd.read_csv("datasets/reporte_1.csv")#,delimiter=';',decimal=',')
df_sismos=pd.read_csv(r'datasets\reporte_LBG.csv')
df_sismos['FECHA - HORA UTC']=df_sismos['Fecha  (UTC)'].astype(str)+' '+df_sismos['Hora  (UTC)'].astype(str)
df_sismos.rename(columns = {'Latitud(°)':'LATITUD (°)', 
                                'Longitud(°)':'LONGITUD (°)',
                                'Profundidad(Km)':'PROF. (Km)',
                                'Magnitud':'MAGNITUD',
                                'Tipo Magnitud':'TIPO MAGNITUD',
                                'Rms(Seg)':'RMS (Seg)',
                                'Gap(°)':'GAP (°)',
                                'Error  Latitud(Km)':'ERROR LATITUD (Km)',
                                'Error  Longitud(Km)':'ERROR LONGITUD (Km)',
                                'Error  Profundidad(Km)':'ERROR PROFUNDIDAD (Km)'}, inplace = True)
df_sismos.drop(['Fecha  (UTC)','Hora  (UTC)'],axis=1,inplace=True)
df_sismos=df_sismos[(df_sismos['PROF. (Km)']<=32)&(df_sismos['PROF. (Km)']>(z_topo.min()*(-1/1000)))&(df_sismos['MAGNITUD']>0)& #Filtros previos
        (df_sismos['LATITUD (°)']>lai)&(df_sismos['LATITUD (°)']<las)
        &(df_sismos['LONGITUD (°)']>loi)&(df_sismos['LONGITUD (°)']<los)]
df_sismos['PROF. (m)']=-df_sismos['PROF. (Km)']*1000 #
df_sismos['ERROR PROFUNDIDAD (m)']=df_sismos['ERROR PROFUNDIDAD (Km)']*1000 #Conversion de km a m
df_sismos['ERROR LONGITUD (°)']=df_sismos['ERROR LONGITUD (Km)']/111.1 #Conversion de km a °
df_sismos['ERROR LATITUD (°)']=df_sismos['ERROR LATITUD (Km)']/111.1 
df_sismos['FECHA - HORA UTC']=df_sismos['FECHA - HORA UTC'].apply(lambda x:pd.to_datetime(x))#Conversion de str a UTCDateTime

#Errores con topografia
df_sismos_err=df_sismos[df_sismos['PROF. (m)']+df_sismos['ERROR PROFUNDIDAD (m)']>(z_topo.min())]
df_sismos_no_err=df_sismos[df_sismos['PROF. (m)']+df_sismos['ERROR PROFUNDIDAD (m)']<=(z_topo.min())]
df_sismos_no_err.loc[:,'ERROR PROFUNDIDAD SUP (m)']=df_sismos_no_err.loc[:,'ERROR PROFUNDIDAD (m)'].copy()
z=[]
for x,y,zhip,zerr in zip(df_sismos_err['LONGITUD (°)'],df_sismos_err['LATITUD (°)'],
                df_sismos_err['PROF. (m)'],df_sismos_err['ERROR PROFUNDIDAD (m)']):
    df_elev=df_topo[(df_topo[0]<(x+0.005))&(df_topo[0]>(x-0.005))&
                (df_topo[1]<(y+0.005))&(df_topo[1]>(y-0.005))]
    d=1
    for x0,y0,z0 in zip(df_elev[0],df_elev[1],df_elev[2]):
        dist=np.sqrt(((x-x0)**2)+((y-y0)**2))
        if dist<d:
            d=dist
            z1=z0
    if z1<=(zhip+zerr):
        z.append(z1-zhip)
    else :
        z.append(zerr)
df_sismos_err.loc[:,'ERROR PROFUNDIDAD SUP (m)']=(np.array(z))
df_sismos=pd.concat([df_sismos_err, df_sismos_no_err])
del df_sismos_err
del df_sismos_no_err

#Inyeccion de H2O
iny=pd.read_csv(r'datasets\inyeccion_geo.csv',delimiter=';',decimal=',')
iny=iny[:-1]


inyec=[]
for name,lon,lat,alt in zip(iny['CAMPO'].apply(lambda x:str(x)),iny['X'],iny['Y'],[2100]*len(iny['CAMPO'])):
    un=dict(
            showarrow=False,
            x=lon,
            y=lat,
            z=alt+10,
            text=name,
            xanchor="left",
            xshift=10,
            opacity=0.7,
            font=dict(
                color="white",
                size=12
            ))
    inyec.append(un)

#PPII
kalei,kicyl1,kibcircles1,kicyl2,kibcircles2,kicyl3,kibcircles3=vol_sus(-73.8566, 7.36551,3902,'Kalé - Investigación','yellow')
kaley,kycyl1,kybcircles1,kycyl2,kybcircles2,kycyl3,kybcircles3=vol_sus(-73.8571014, 7.3647799,2618.232,'Kalé - Inyector','blue')
plai,picyl1,pibcircles1,picyl2,pibcircles2,picyl3,pibcircles3=vol_sus(-73.89389, 7.2572498,3227.8,'Platero - Investigación','red')
play,pycyl1,pybcircles1,pycyl2,pybcircles2,pycyl3,pybcircles3=vol_sus(-73.8944016, 7.25667,2325.6,'Platero - Inyector','green')

kalec,_,_,_,_,_,_=vol_sus(-73.8570023, 7.3647499,2325.6,'Kalé - Captador','orange')
plac,_,_,_,_,_,_=vol_sus(-73.8943024, 7.2566800,2325.6,'Platero - Captador','gold')

#Estaciones sismologicas
df_sta_vmm=pd.read_csv('datasets/VMM_STA.csv',delimiter=';',decimal=',')
df_sta_lom=pd.read_csv('datasets//LOMA_STA.csv',delimiter=';',decimal=',')

STA_VMM = go.Scatter3d(
    x=df_sta_vmm['LONGITUD'],
    y=df_sta_vmm['LATITUD'],
    z=df_sta_vmm['ALTITUD (msnm)']+10, #Para sobresalir de la topografía
    mode='markers',
    marker_symbol='diamond',
    name="Estación sismologica VMM",
    hovertemplate ='Longitud:'+df_sta_vmm['LONGITUD'].apply(lambda x:str(x))+'°'+'<br>'+
                    'Latitud:'+df_sta_vmm['LATITUD'].apply(lambda x:str(x))+'°'+'<br>'+
                    'Elevacion:'+df_sta_vmm['ALTITUD (msnm)'].apply(lambda x:str(x))+'msnm <br>'+
                    'Nombre de la estación:'+df_sta_vmm['NOMBRE ESTACIÓN']+'°'+'<br>'+
                    'Código:'+df_sta_vmm['CODIGO']+'<br>'+
                    'Agencia:'+df_sta_vmm['AGENCIA']+'<br>'+
                    'Fecha de instalación:'+df_sta_vmm['FECHA DE INSTALACIÓN'].apply(lambda x:str(x))+'<br>'+
                    'Fecha de retiro:'+df_sta_vmm['FECHA DE RETIRO'].apply(lambda x:str(x))+'<br>'+
                    'Estado:'+df_sta_vmm['ESTADO'],
    marker=dict(
        size=4,
        color='blueviolet'
    ),
    showlegend=False
)
STA_LOM = go.Scatter3d(
    x=df_sta_lom['LONGITUD'],
    y=df_sta_lom['LATITUD'],
    z=df_sta_lom['ALTITUD (msnm)']+10,
    mode='markers',
    marker_symbol='diamond',
    name="Estación sismologica la Loma, Cesar",
    hovertemplate ='Longitud:'+df_sta_lom['LONGITUD'].apply(lambda x:str(x))+'°'+'<br>'+
                    'Latitud:'+df_sta_lom['LATITUD'].apply(lambda x:str(x))+'°'+'<br>'+
                    'Elevacion:'+df_sta_lom['ALTITUD (msnm)'].apply(lambda x:str(x))+'msnm <br>'+
                    'Nombre de la estación:'+df_sta_lom['NOMBRE ESTACIÓN']+'°'+'<br>'+
                    'Código:'+df_sta_lom['CODIGO']+'<br>'+
                    'Agencia:'+df_sta_lom['AGENCIA']+'<br>'+
                    'Fecha de instalación:'+df_sta_lom['FECHA DE INSTALACIÓN'].apply(lambda x:str(x))+'<br>'+
                    'Fecha de retiro:'+df_sta_lom['FECHA DE RETIRO'].apply(lambda x:str(x))+'<br>'+
                    'Estado:'+df_sta_lom['ESTADO'],
    marker=dict(
        size=6,
        color='blueviolet'
    ),
    showlegend=False
)

#Rios
df_rivers=pd.read_csv('datasets\drenajes_SIM.csv')

#Cargar datos de pozos
df_pozos=pd.read_csv('datasets/pozos.csv',usecols=['lon', 'lat', 'UWI', 'WELL_NAME', 
'DEPARTAMEN', 'WELL_COU_1', 'WELL_TVD', 'WELL_KB_EL',
       'ROTARY_ELE', 'WELL_DRILL', 'WELL_GROUN', 'FIELD_ABRE',
       'CONTRATO', 'WELL_SPUD_', 'COORD_QUAL', 'COMMENT_', 'WELL_COMPL', 'WELL_STA_1',
       'WELLTYPE', 'FECHA_ACTU',
       'OPERATOR_W', 'COMPANY_CO', 'z'])
Pozos = go.Scatter3d(
    x=df_pozos['lon'],
    y=df_pozos['lat'],
    z=df_pozos['z']+15,
    mode='markers',
    name="Pozo petrolífero",
    hovertemplate ='Longitud:'+df_pozos['lon'].apply(lambda x:str(x))+'°'+'<br>'+
                    'Latitud:'+df_pozos['lat'].apply(lambda x:str(x))+'°'+'<br>'+
                    'Elevacion:'+df_pozos['z'].apply(lambda x:str(x))+'msnm <br>'+
                    'UWI:'+df_pozos['UWI']+'°'+'<br>'+
                    'Nombre del pozo:'+df_pozos['WELL_NAME']+'<br>'+
                    'Departamento:'+df_pozos['DEPARTAMEN']+'<br>'+
                    'Municipio:'+df_pozos['WELL_COU_1']+'<br>'+
                    'Tipo:'+df_pozos['WELLTYPE']+'<br>'+
                    'Operador:'+df_pozos['OPERATOR_W']+'<br>'+
                    'Compañia:'+df_pozos['COMPANY_CO'],
    marker=dict(
        size=1.5,
        color='black'
    )
)

#Rezumaderos
df_rezumaderos=pd.read_csv('datasets\REZUMADEROS_WGS84_SIM.txt',decimal=',',delimiter=';')
rez_txt=('Longitud:'+df_rezumaderos['X'].apply(lambda x:str(x))+'°'+
            '<br>'+'Latitud:'+df_rezumaderos['Y'].apply(lambda x:str(x))+'°'+
            '<br>'+'Elevacion:'+df_rezumaderos['Z'].apply(lambda x:str(x))+'msnm'+
            '<br>'+'Tipo:'+df_rezumaderos['TIPO'].apply(lambda x:str(x))+
            '<br>'+'Autor:'+df_rezumaderos['AUTOR'].apply(lambda x:str(x))+
            '<br>'+'Empresa:'+df_rezumaderos['EMPRESA'].apply(lambda x:str(x))+
            '<br>'+'Formación:'+df_rezumaderos['FORMACION_'].apply(lambda x:str(x))+
            '<br>'+'Tipo secundario:'+df_rezumaderos['TIPO_2'].apply(lambda x:str(x))+
            '<br>'+'Departamento:'+df_rezumaderos['DPTOS_NOMB'].apply(lambda x:str(x))+
            '<br>'+'Capital:'+df_rezumaderos['CAPITAL'].apply(lambda x:str(x)))
rez = go.Scatter3d(
    x=df_rezumaderos['X'],
    y=df_rezumaderos['Y'],
    z=df_rezumaderos['Z'],
    mode='markers',
    name="Rezumaderos",
    hovertemplate =rez_txt,
    marker=dict(
        size=5,
        color='gray'
    ),
    showlegend=False
)
#Poblaciones
df_poblaciones=pd.read_csv('datasets/poblaciones.csv',usecols=['Name','lon','lat','outputSRTM1'])
Poblaciones = go.Scatter3d(
    x=df_poblaciones['lon'],
    y=df_poblaciones['lat'],
    z=df_poblaciones['outputSRTM1']+10, #Conseguir alturas
    mode='markers',
    name="Población",
    marker_symbol='square',
    hovertemplate =df_poblaciones['Name'],
    marker=dict(
        size=6,
        color='red'
    ),
    textposition="bottom right"
)
Pobl=[]
for name,lon,lat,alt in zip(df_poblaciones['Name'],df_poblaciones['lon'],
df_poblaciones['lat'],df_poblaciones['outputSRTM1']):
    un=dict(
            showarrow=False,
            x=lon,
            y=lat,
            z=alt+10,
            text=name,
            xanchor="left",
            xshift=10,
            opacity=0.7,
            font=dict(
                color="white",
                size=12
            ))
    Pobl.append(un)

#Carreteras
roads=pd.read_csv('datasets\Via_WGS84_SIM.txt',delimiter=';',decimal=',')
roads  =roads[(roads['LATITUD']>lai)&(roads['LATITUD']<las)&(roads['LONGITUD']>loi)&(roads['LONGITUD']<los)]

# ls_x,ls_y,ls_z=lin_list('datasets/fallas.csv') #Fallas
fallas=pd.read_csv('datasets/fallas_SIM.csv',decimal=',')
fallas['X']=fallas['X'].apply(lambda x:float(x))
fallas['Y']=fallas['Y'].apply(lambda x:float(x))
fallas['Z']=fallas['Z'].apply(lambda x:float(x))
fallas_1=pd.read_csv('datasets/fallas_1_SIM.csv')
fallas_1=fallas_1.drop_duplicates(subset=['LINE_ID'])

# ls_x_f,ls_y_f,ls_z_f=lin_list('datasets/campos.csv') #Campos
campet=pd.read_csv('datasets\campos_SIM.csv',decimal=',')
campet['X']=campet['X'].apply(lambda x:float(x))
campet['Y']=campet['Y'].apply(lambda x:float(x))
campet['Z']=campet['Z'].apply(lambda x:float(x))
campet_1=pd.read_csv('datasets/campos_1_SIM.csv')
campet_1=campet_1.drop_duplicates(subset=['LINE_ID'])

# ls_x_s,ls_y_s,ls_z_s=lin_list('datasets/lineas.csv') #Lineas sismicas
linsis=pd.read_csv('datasets\lineas_SIM.csv',decimal=',')
linsis['X']=linsis['X'].apply(lambda x:float(x))
linsis['Y']=linsis['Y'].apply(lambda x:float(x))
linsis['Z']=linsis['Z'].apply(lambda x:float(x))
linsis_1=pd.read_csv('datasets/lineas_1_SIM.csv')
linsis_1=linsis_1.drop_duplicates(subset=['LINE_ID'])



Eoceno=geology('datasets/DISCORDANCIA_EOCENO.txt','pink','red','Discordancia del Eoceno Medio')
Colorado=geology('datasets/TOPE_COLORADO.txt','darkblue','aquamarine','Tope Formación Colorado')
Mugrosa=geology('datasets/TOPE_MUGROSA.txt','green','greenyellow','Tope Formación Mugrosa')
Chorros=geology('datasets/TOPE_CHORROS.txt','orangered','yellow','Tope Grupo Chorros')
Real=geology('datasets/BASE_CUATERNARIO.txt','purple','pink','Tope Grupo Real')





df_new=pd.read_csv('datasets/UN_CRN_COLORS.csv',index_col=None)


SISMICA=img_3d("ANH-TR-2006-04-A","assets\perfil_2_sintexto.jpg",-74.115,7.58,-72.954,6.806,4300,-20000)

#Texto de imagenes
df_andina=pd.read_csv('datasets\Coordenadas_textos_perfil_trasandina.csv',delimiter=';')
txts_p=[]
for name,lon,lat,alt in zip(df_andina['texto'],df_andina['x'],df_andina['y'], df_andina['z']):
    un=dict(
            showarrow=False,
            x=lon,
            y=lat,
            z=alt,
            text=name,
            xshift=0,
            opacity=0.7,
            font=dict(
                color='white',
                size=12
            ))
    txts_p.append(un)

#exp = open("datasets\Explicacion_modelo3d.txt", "r")
exp=io.open("datasets\Explicacion_modelo3d.txt", mode="r", encoding="utf-8")
ls_k=[html.H2('¿Cómo funciona el modelo tridimensional del Valle Medio del Magdalena?', className="card-text")]
for i in exp:
    ls_k.append(html.H6(i, className="card-text"))

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.SUPERHERO])
app.config['suppress_callback_exceptions'] = True
#Cargars los datos

card_main=dbc.Card(
    dbc.CardBody(
        
            [dbc.Nav([
                dbc.NavLink("Inicio", href="https://www.centrodetransparenciappii.org/", active="exact"),
                dbc.NavLink("Modelo 3D VMM", href="", active="exact"),
                dbc.NavLink("Semáforo sísmico", href="https://pinguinodigital.com/wp-content/uploads/2020/08/pagina-en-construcci%C3%B3n1.jpg", active="exact"),
            ]),
            html.H2("Modelo Tridimensional de Sismicidad en el Valle Medio del Magdalena", className="card-title"),
            html.H4("Transparencia:", className="card-subtitle"),
            html.H5("Topografia:", className="card-subtitle"),
            dcc.Slider(
                id='TOPO',
                min=0,
                max=1,
                step=0.1,
                value=1,
                tooltip={"placement": "bottom", "always_visible": False}),
                #Condicionales de geologia-Tope Grupo real
                html.Div(id='GREAL', children=[
                        html.H5("Tope Grupo Real:", className="card-subtitle"),
                        # Create element to hide/show, in this case a slider
                        dcc.Slider(id='TGREAL',
                                min=0,
                                max=1,
                                step=0.1,
                                value=1,
                                tooltip={"placement": "bottom", "always_visible": False})

                    ], style= {'display': 'none'}),
                #Condicionales de geologia-Tope Grupo real
                html.Div(id='COLORADO', children=[
                        html.H5("Tope Formación Colorado:", className="card-subtitle"),
                        # Create element to hide/show, in this case a slider
                        dcc.Slider(id='TCOLORADO',
                                min=0,
                                max=1,
                                step=0.1,
                                value=1,
                                tooltip={"placement": "bottom", "always_visible": False})

                    ], style= {'display': 'none'}),
                #Condicionales de geologia-Tope Grupo real
                html.Div(id='MUGROSA', children=[
                        html.H5("Tope Formación Mugrosa:", className="card-subtitle"),
                        # Create element to hide/show, in this case a slider
                        dcc.Slider(id='TMUGROSA',
                                min=0,
                                max=1,
                                step=0.1,
                                value=1,
                                tooltip={"placement": "bottom", "always_visible": False})

                    ], style= {'display': 'none'}),
                #Condicionales de geologia-Tope Grupo real
                html.Div(id='CHORROS', children=[
                        html.H5("Tope Formación Chorros:", className="card-subtitle"),
                        # Create element to hide/show, in this case a slider
                        dcc.Slider(id='TCHORROS',
                                min=0,
                                max=1,
                                step=0.1,
                                value=1,
                                tooltip={"placement": "bottom", "always_visible": False})

                    ], style= {'display': 'none'}),
                #Condicionales de geologia-Tope Grupo real
                html.Div(id='EOCMED', children=[
                        html.H5("Discordancia del Eoceno Medio:", className="card-subtitle"),
                        # Create element to hide/show, in this case a slider
                        dcc.Slider(id='TEOCMED',
                                min=0,
                                max=1,
                                step=0.1,
                                value=1,
                                tooltip={"placement": "bottom", "always_visible": False})

                    ], style= {'display': 'none'}),
           #Fin condicionales
            html.H4("Exageración vertical:", className="card-subtitle"),
            dcc.Slider(
                id='EXG',
                min=1,
                max=10,
                step=1,
                value=2,
                tooltip={"placement": "bottom", "always_visible": False}),
            html.H4("Magnitudes:", className="card-subtitle"),
                    dcc.RangeSlider(
                id='MAGN',
                min=df_sismos['MAGNITUD'].min(),
                max=df_sismos['MAGNITUD'].max(),
                step=0.1,
                value=[df_sismos['MAGNITUD'].min(), df_sismos['MAGNITUD'].max()],
                # marks=None,
                allowCross=False,
                tooltip={"placement": "bottom", "always_visible": False}
            ),
            html.H4("Profundidad (m):", className="card-subtitle"),
            dcc.RangeSlider(
                id='DEPTH',
                min=df_sismos['PROF. (m)'].min(),
                max=df_sismos['PROF. (m)'].max(),
                step=100,
                value=[df_sismos['PROF. (m)'].min(), df_sismos['PROF. (m)'].max()],
                # marks=None,
                allowCross=False,
                tooltip={"placement": "bottom", "always_visible": False}
            ),
            html.H4("Fecha:", className="card-subtitle"),
            dcc.DatePickerRange(
                    id='DATE',
                    start_date_placeholder_text="Start Date",
                    end_date_placeholder_text="End Date",
                    calendar_orientation='horizontal',
                    start_date=df_sismos['FECHA - HORA UTC'].min(),
                    end_date=df_sismos['FECHA - HORA UTC'].max(),
                    day_size=30,
                    min_date_allowed=df_sismos['FECHA - HORA UTC'].min(),
                    max_date_allowed=df_sismos['FECHA - HORA UTC'].max(),
                    persistence=True,
                    #initial_visible_month=df_sismos['FECHA - HORA UTC'].min(),
                    reopen_calendar_on_clear=False
                ),
                html.H4("Perfiles:", className="card-subtitle"),
                html.H5("Punto 1 (Longitud-Latitud)", className="card-subtitle"),
                dcc.Input(id="Longitud 1", type="number", placeholder="Longitud 1", min=loi, max=los, step=0.01,style={'marginRight':'10px'},value=loi),
                dcc.Input(id="Latitud 1", type="number", placeholder="Latitud 1", min=lai, max=las, step=0.01, debounce=True,value=lai),
                html.H5("Punto 2 (Longitud-Latitud)", className="card-subtitle"),
                dcc.Input(id="Longitud 2", type="number", placeholder="Longitud 2", min=loi, max=los, step=0.01,value=los,style={'marginRight':'10px'}),
                dcc.Input(id="Latitud 2", type="number", placeholder="Latitud 2", min=lai, max=las, step=0.01,value=las,debounce=True),
            html.H4("Variables de sismicidad desplegadas:", className="card-subtitle"),
            dcc.Dropdown(id='SEISMO',
                        placeholder="Variables a desplegar...",
                        style={'color': 'black'},
                        options=[
                            {'label': 'Localización', 'value': 'LOC'},
                            {'label': 'Fecha', 'value': 'FEC'},
                            {'label': 'Magnitud', 'value': 'MAG'},
                            {'label': 'RMS', 'value': 'RMS'},
                            {'label': 'Errores', 'value': 'ERR'},
                            {'label': 'Mostrar sismicidad', 'value': 'SISM'}
                        ],
                        value=['LOC', 'FEC','MAG','RMS','ERR','SISM'],
                        multi=True
                    ),
            html.H4("________________________________________", className="card-subtitle"),
            html.H4("PPII:", className="card-subtitle"),
            dcc.Dropdown(id='PPII',
                        placeholder="Variables a desplegar...",
                        style={'color': 'black'},
                        options=[
                            {'label': ' Pozo Kalé - Investigación (ANH)', 'value': 'KALEi'},
                            {'label': ' Pozo Kalé - Inyector (ANH)', 'value': 'KALEy'},
                            {'label': ' Pozo Kalé - Captador (ANH)', 'value': 'KALEc'},
                            {'label': ' Pozo Platero - Investigación (ANH)', 'value': 'PLATEROi'},
                            {'label': ' Pozo Platero - Inyector (ANH)', 'value': 'PLATEROy'},
                            {'label': ' Pozo Platero - Captador (ANH)', 'value': 'PLATEROc'},
                            
                        ],
                        value=[],
                        multi=True
                    ),
            html.H4("________________________________________", className="card-subtitle"),
            html.H4("Cartografía base:", className="card-subtitle"),
            dcc.Dropdown(id='CART',
                        placeholder="Variables a desplegar...",
                        style={'color': 'black'},
                        options=[
                            {'label': ' Barras de Error (SGC)', 'value': 'ERROR'},
                            {'label': ' Estaciones sismológicas (SGC)', 'value': 'STA'},
                            {'label': ' Poblaciones (UNAL-ANH-MINCIENCIAS)', 'value': 'POB'},
                            {'label': ' Drenajes (IGAC)', 'value': 'RIV'},
                            {'label': ' Vias (IGAC)', 'value': 'VIA'},
                            {'label': ' Perfil', 'value': 'PER'},
                            
                        ],
                        value=[],
                        multi=True
                    ),
            html.H4("________________________________________", className="card-subtitle"),
            html.H4("Información complementaria:", className="card-subtitle"),
            dcc.Dropdown(id='PETRO',
                        placeholder="Variables a desplegar...",
                        style={'color': 'black'},
                        options=[
                            {'label': ' Pozos petrolíferos (UNAL-ANH-MINCIENCIAS)', 'value': 'POZO'},
                            {'label': ' Campos petrolíferos (UNAL-ANH-MINCIENCIAS)', 'value': 'FIELD'},
                            {'label': ' Trazo en superficie de líneas sísmicas (UNAL-ANH-MINCIENCIAS)', 'value': 'LIN'},
                            {'label': ' Rezumaderos (ANH)', 'value': 'REZ'},
                            {'label': ' ANH-TR-2006-04-A (ANH)', 'value': 'SEIS'},
                            {'label': ' Inyección de agua', 'value': 'H2O'}
                            
                        ],
                        value=[],
                        multi=True
                    ),
#------------Condicional campos--------
                html.Div(id='INY', children=[
                        html.H5("Campos (inyección de agua para recobro mejorado):", className="card-subtitle"),
                        # Create element to hide/show, in this case a slider
                        dcc.Dropdown(id='TINY',
                                    placeholder="Campo",
                                    style={'color': 'black'},
                                    options=[{'label':x,'value':x} for x in iny['CAMPO']],
                                    value='LA CIRA',
                                    multi=False,
                                    clearable=False
                                )

                    ], style= {'display': 'none'}),
#----------------------------------
            html.H4("________________________________________", className="card-subtitle"),
            html.H4("Geología:", className="card-subtitle"),
            dcc.Dropdown(id='GEOL',
                        placeholder="Variables a desplegar...",
                        style={'color': 'black'},
                        options=[
                            {'label': ' Fallas Geológicas (SGC)', 'value': 'FALL'},
                            {'label': ' Tope Grupo Real (UNAL-ANH-MINCIENCIAS)', 'value': 'REAL'},
                            {'label': ' Tope Formación Colorado (UNAL-ANH-MINCIENCIAS)', 'value': 'COL'},
                            {'label': ' Tope Formación Mugrosa (UNAL-ANH-MINCIENCIAS)', 'value': 'MUG'},
                            {'label': ' Tope Grupo Chorros (UNAL-ANH-MINCIENCIAS)', 'value': 'CHO'},
                            {'label': ' Discordancia del Eoceno Medio (UNAL-ANH-MINCIENCIAS)', 'value': 'EOC'},
                            {'label': ' Geología superficial (SGC)', 'value': 'GEO'},
                            
                        ],
                        value=[],
                        multi=True
                    ),
            dcc.Markdown('''
                    * ANH: Agencia Nacional de Hidrocarburos
                    * MINCIENCIAS: Ministerio de Ciencia Tecnología e Innovación
                    * SGC: Servicio Geológico Colombiano
                    * UNAL: Universidad Nacional de Colombia
                '''),
                     dbc.CardImg(src="assets\logos.png", bottom=True, alt='Logos_convenio_tripartito',)    
                 ,], 
        ),
    color="secondary",   # https://bootswatch.com/default/ for more card colors
    inverse=True,   # change color of text (black or white)
    # outline=False,  # True = remove the block colors from the background and header,
    style={"overflow": "scroll","width": "28rem",'height':'40rem'},
)

card_graph = dbc.Card(
        dcc.Graph(id='3d_model', figure={}), body=True,color="dark",style={'height':'40rem'}
)

card_graph_profile = dbc.Card(
        dcc.Graph(id='Model_profile', figure={}), body=True,color="dark",
)

card_iny_graph = dbc.Card(
        dcc.Graph(id='Iny_graph', figure={}), body=True,color="dark",
)

references=[
        # html.H2("Explicacion del Modelo 3d", className="card-title"),
        # html.H6(exp, className="card-text"),
        html.H2("Referencias", className="card-title"),
        html.H6("Agencia Nacional de Hidrocarburos - ANH & Servicio Geológico Colombiano - SGC (2016). Informe final del Convenio interadministrativo 194 ANH-014 SGC, entre la Agencia Nacional de Hidrocarburos y el Servicio Geológico Colombiano.", 
            className="card-text"),
        html.H6("Agencia Nacional de Hidrocarburos - ANH (2010). Mapa de Rezumaderos. Información Geológica y Geofísica. https://www.anh.gov.co/Informacion-Geologica-y-Geofisica/Estudios-Integrados-y-Modelamientos/Paginas/MAPA-DE-REZUMADEROS.aspx", 
            className="card-text"),         
        html.H6("Ángel-Martínez, C.E., Prieto-Gómez, G.A., Cristancho-Mejía, F., Sarmiento-Orjuela, A.M., Vargas-Quintero, J.A., Delgado-Mateus, C.J., Torres-Rojas, E., Castelblanco-Ossa, C.A., Camargo-Rache, G.L., Amazo-Gómez, D.F., Cipagauta-Mora, J.B., Lucuara-Reyes, E.D., Ávila-López, K.L. Fracica-González, L.R., Martín-Ravelo, A.S., Atuesta-Ortiz, D.A., Gracía-Romero, D.F., Triviño Cediel , R.J., Jaimes Villarreal, V.N., y Alarcón Rodríguez, W.F.(2021). Proyecto MEGIA: Modelo Geológico-Geofísico del Valle Medio del Magdalena. Producto No. 5. Bogotá: 192 pp.", 
            className="card-text"),
        html.H6("Dionicio, V., Mercado, O. y Lizarazo, M. (2020). Semáforo para el monitoreo sísmico durante el desarrollo de los proyectos piloto de investigación integral en yacimientos no convencionales de hidrocarburos en Colombia. Bogotá: Servicio Geológico Colombiano.", 
            className="card-text"),
        html.H6("Gómez, J. & Montes, N.E., compiladores. 2020. Mapa Geológico de Colombia 2020. Escala 1:1 000 000. Servicio Geológico Colombiano, 2 hojas. Bogotá.​", 
            className="card-text"),
        html.H6("Instituto Geográfico Agustin Codazzi - IGAC (2019). Base de datos vectorial básica. Colombia. Escala 1:100.000. Colombia en Mapas. https://www.colombiaenmapas.gov.co/#", 
            className="card-text"),
        html.H6("Servicio Geológico Colombiano. (2021). Banco de Información Petrolera. https://srvags.sgc.gov.co/JSViewer/GEOVISOR_BIP/", 
            className="card-text"),
        html.H6("Servicio Geológico Colombiano. (2022). Catálogo línea base de sismicidad: Valle Medio del Magdalena y La Loma Cesar. http://bdrsnc.sgc.gov.co/paginas1/catalogo/Consulta_Valle_Medio/valle_medio.php", 
            className="card-text"),
         html.H6("Tozer, B, Sandwell, D. T., Smith, W. H. F., Olson, C., Beale, J. R., & Wessel, P. (2019). Global bathymetry and topography at 15 arc sec: SRTM15+. Distributed by OpenTopography. https://doi.org/10.5069/G92R3PT9. Accessed: 2022-02-10", 
            className="card-text"),
    ]

ls_k.extend(references)

card_references=dbc.Card(
    dbc.CardBody(ls_k
    ))

app.layout = html.Div([
    dbc.Row([dbc.Col(card_main, width=4),
             dbc.Col(card_graph, width=8),
             dbc.Col(card_graph_profile, width=12),
             dbc.Col(card_iny_graph, width=12),
             dbc.Col(card_references, width=12)], 
             justify="start"), 
             
])


@app.callback(
     [dash.dependencies.Output(component_id='3d_model', component_property='figure'),
      dash.dependencies.Output(component_id='DATE', component_property='initial_visible_month'),
      dash.dependencies.Output(component_id='GREAL', component_property='style'),
      dash.dependencies.Output(component_id='COLORADO', component_property='style'),
      dash.dependencies.Output(component_id='MUGROSA', component_property='style'),
      dash.dependencies.Output(component_id='CHORROS', component_property='style'),
      dash.dependencies.Output(component_id='EOCMED', component_property='style'),
      dash.dependencies.Output(component_id='INY', component_property='style')
      ],



    [dash.dependencies.Input(component_id='TOPO', component_property='value'),
     dash.dependencies.Input(component_id='EXG', component_property='value'),
     dash.dependencies.Input(component_id='DATE', component_property='start_date'),
     dash.dependencies.Input(component_id='DATE', component_property='end_date'),
     dash.dependencies.Input(component_id='MAGN', component_property='value'),
     dash.dependencies.Input(component_id='DEPTH', component_property='value'),
     dash.dependencies.Input(component_id='SEISMO', component_property='value'),
     dash.dependencies.Input(component_id='PPII', component_property='value'),
     dash.dependencies.Input(component_id='CART', component_property='value'),
     dash.dependencies.Input(component_id='PETRO', component_property='value'),
     dash.dependencies.Input(component_id='GEOL', component_property='value'),
     dash.dependencies.Input(component_id='Longitud 1', component_property='value'),
     dash.dependencies.Input(component_id='Longitud 2', component_property='value'),
     dash.dependencies.Input(component_id='Latitud 1', component_property='value'),
     dash.dependencies.Input(component_id='Latitud 2', component_property='value'),
     dash.dependencies.Input(component_id='TGREAL', component_property='value'),
     dash.dependencies.Input(component_id='TCOLORADO', component_property='value'),
     dash.dependencies.Input(component_id='TMUGROSA', component_property='value'),
     dash.dependencies.Input(component_id='TCHORROS', component_property='value'),
     dash.dependencies.Input(component_id='TEOCMED', component_property='value'),
     ])

def update_figure(TOPO,EXG,START_DATE,END_DATE,MAGN,DEPTH,SEISMO,PPII,CART,PETRO,GEOL,x0,x1,y0,y1,
                        TGREAL,TCOLORADO,TMUGROSA,TCHORROS,TEOCMED):
        sub=None
        fig=go.Figure()
        
        if np.isin('H2O', PETRO):
            INYO={'display': 'block'}
            iny=pd.read_csv(r'datasets\inyeccion_geo.csv',delimiter=';',decimal=',')
            iny=iny[:-1]
            ls_form=[False]*len(iny['CAMPO'])
            ls_form[-1]=True
            # try:
            for i,cond in zip(iny['CAMPO'],ls_form):
                    inyc=iny[iny['CAMPO']==i]
                    fig.add_trace(go.Scatter3d(x=[float(inyc['X'])]*2, y=[float(inyc['Y'])]*2, z=[0,-1*float(inyc['prof'])],
                                hovertemplate=[inyc['CAMPO'].apply(lambda x:str(x))+'<br>'
                                        'Pozos:'+inyc['POZOS'].apply(lambda x:str(x))+'<br>'
                                        'BBL:'+inyc['TOTAL_bbl'].apply(lambda x:str(x))]*2,mode='lines',name=str(i),
                                        line=dict(color=inyc['TOTAL_bbl'],
                                    width=20,colorscale='Jet',
                                cmax=((iny['TOTAL_bbl'])).max(),
                                cmin=((iny['TOTAL_bbl'])).min(),
                                #showscale=cond,
                                #colorbar={"title": 'Volumen de <br>inyección (BBL)','x': 1.25}
                                )
                                ,showlegend=False),)
        else:
            INYO={'display': 'none'}
        if np.isin('GEO', GEOL):
            if TOPO==0:
                DISM=0
            else:
                DISM=0.01
            topog.colorscale=['black','black']
            topog.opacity=TOPO-DISM
            fig.add_trace(topog,row=sub,col=sub)
            directory = 'datasets\CSV_UNIDADES'
            for filename in os.scandir(directory):
                if filename.is_file():
                    name=(str(filename.path).split('\\'))[-1]
                    name=name.replace('.csv','')
                    name=name.replace('_','?')
                    df_1=df_new[df_new['SimboloUC']==name]
                    text='Edad: '+np.array(df_1['Edad'].apply(lambda x:str(x)))[0]+'<br>Descripción: '+np.array(df_1['Descripcio'].apply(lambda x:str(x)))[0]#+'<br>UGIntegrad: '+np.array(df_1['UGIntegrad'].apply(lambda x:str(x)))[0]
                    text=str(text)
                    fig.add_trace(geology_super(filename.path,np.array(df_1['Color'].apply(lambda x:str(x)))[0],name,text,TOPO))
            
        else:
            topog.colorscale=['green','greenyellow','saddlebrown','saddlebrown','saddlebrown','saddlebrown','snow','snow']
            topog.opacity=TOPO
            fig.add_trace(topog,row=sub,col=sub)
        df_sismos_1=df_sismos[(df_sismos['FECHA - HORA UTC']<=END_DATE)&(df_sismos['FECHA - HORA UTC']>=START_DATE)&
        (df_sismos['MAGNITUD']>=MAGN[0])&(df_sismos['MAGNITUD']<=MAGN[1])
        &(df_sismos['PROF. (m)']>=DEPTH[0])&(df_sismos['PROF. (m)']<=DEPTH[1])]
        text=text_scatter(SEISMO,df_sismos_1)
        if np.isin('SISM',SEISMO):
            vis=True
        else :
            vis=False
        if np.isin('ERROR', CART):
            err=True
        else:
            err=False
        fig.add_trace(go.Scatter3d(
            x=df_sismos_1['LONGITUD (°)'],y=df_sismos_1['LATITUD (°)'],z=df_sismos_1['PROF. (m)'],mode='markers',
            
            marker=dict(
                size=(df_sismos_1['MAGNITUD'])**2,
                color=df_sismos_1['PROF. (m)'],                # set color to an array/list of desired values
                colorscale='Jet',   # choose a colorscale
                opacity=0.8,
                cmax=df_sismos['PROF. (m)'].max(),
                cmin=-32000,
                #showscale=True,
                #colorbar={"title": 'Profundidad del <br> sismo (m)',
                #    "orientation": 'h'},
            ),
            error_x=dict(
                array=df_sismos_1['ERROR LONGITUD (°)'],                # set color to an array/list of desired values
                color='red',   # choose a colorscale
                symmetric=True,
                width=0.01,
                visible=err
            ),
            error_y=dict(
                array=df_sismos_1['ERROR LATITUD (°)'],                # set color to an array/list of desired values
                color='red',   # choose a colorscale
                symmetric=True,
                width=0.01,
                visible=err
            ),
            error_z=dict(
                array=df_sismos_1['ERROR PROFUNDIDAD SUP (m)'], 
                arrayminus=df_sismos_1['ERROR PROFUNDIDAD (m)'] ,             
                color='red',   # choose a colorscale
                symmetric=False,
                width=0.01,
                visible=err
            ),
            hovertemplate=text,
                name='Sismos',
                showlegend=False,
                visible=vis
                ),row=sub,col=sub)
        # if np.isin('KALE', CART):
        #     fig.add_trace(kale)
        if np.isin('KALEi', PPII):
            fig.add_traces([kalei,kicyl1,kibcircles1,kicyl2,kibcircles2,kicyl3,kibcircles3])
        if np.isin('KALEy', PPII):
            fig.add_traces([kaley,kycyl1,kybcircles1,kycyl2,kybcircles2,kycyl3,kybcircles3])
        if np.isin('KALEc', PPII):
            fig.add_trace(kalec)
        if np.isin('PLATEROi', PPII):
            fig.add_traces([plai,picyl1,pibcircles1,picyl2,pibcircles2,picyl3,pibcircles3])
        if np.isin('PLATEROy', PPII):
            fig.add_traces([play,pycyl1,pybcircles1,pycyl2,pybcircles2,pycyl3,pybcircles3])
        if np.isin('PLATEROc', PPII):
            fig.add_trace(plac)

        if np.isin('RIV', CART):
            for i in df_rivers['DRENAJE'].unique():
                riv=df_rivers[df_rivers['DRENAJE']==i]
                fig.add_trace(go.Scatter3d(z=riv['Z'], x=riv['X'], y=riv['Y'],mode='markers',
                name=str(i),marker_symbol='square',marker=dict(color='aqua',size=3),showlegend=False,hovertemplate=str(i)))
        if np.isin('STA', CART):
            fig.add_trace(STA_VMM)
            fig.add_trace(STA_LOM)
        if np.isin('VIA', CART):
            for i in roads['GLOBALID'].unique():
                f=roads[roads['GLOBALID']==i]
                fig.add_trace(go.Scatter3d(x=f['LONGITUD'], y=f['LATITUD'], z=f['ELEVACION'],
                hovertemplate=str(i),mode='lines',name='Via',line=dict(color='yellow',width=2),showlegend=False),)
        if np.isin('POZO', PETRO):
            fig.add_trace(Pozos)

        if np.isin('REZ', PETRO):
            fig.add_trace(rez)
        if np.isin('POB', CART):
                fig.add_trace(Poblaciones)
                fig.update_layout(
                scene=dict(
                annotations=Pobl),
                overwrite=False)
        if np.isin('FIELD', PETRO):
            for i in campet['LINE_ID'].unique():
                f=campet[campet['LINE_ID']==i]
                attr=campet_1[campet_1['LINE_ID']==i]
                nom='Compañia:'+np.array(attr['Compañia'])[0]+'<br>Estado:'+np.array(attr['Estado'])[0]+'<br>Información:'+str(np.array(attr['INFO'])[0])
                try:
                    tip='Campo petrolífero '+np.array(attr['Campo'])[0]
                except:
                    tip='_'
                fig.add_trace(go.Scatter3d(x=f['X'], y=f['Y'], z=f['Z'],
                                hovertemplate=nom,
                                mode='lines',
                                name=tip,line=dict(color='black',width=3),showlegend=False),)
                fig.update_layout(
                        scene=dict(
                        annotations=inyec),
                overwrite=False)
        if np.isin('LIN', PETRO):
            for i in linsis['LINE_ID'].unique():
                f=linsis[linsis['LINE_ID']==i]
                attr=linsis_1[linsis_1['LINE_ID']==i]
                try:
                    nom='ssGmStNm:'+np.array(attr['ssGmStNm'])[0]+'<br>Proyecto:'+np.array(attr['project'])[0]
                except:
                    nom='_'
                try:
                    tip=np.array(attr['owtype'])[0]
                except:
                    tip='_'
                fig.add_trace(go.Scatter3d(x=f['X'], y=f['Y'], z=f['Z'],
                                hovertemplate=nom,
                                mode='lines',
                                name=tip,line=dict(color='blue',width=3),showlegend=False),)
        if np.isin('FALL', GEOL):
            for i in fallas['LINE_ID'].unique():
                f=fallas[fallas['LINE_ID']==i]
                attr=fallas_1[fallas_1['LINE_ID']==i]
                try:
                    nom=np.array(attr['NombreFall'])[0]
                except:
                    nom='_'
                try:
                    tip=np.array(attr['Tipo'])[0]
                except:
                    tip='_'
                fig.add_trace(go.Scatter3d(x=f['X'], y=f['Y'], z=f['Z'],
                                hovertemplate=nom,
                                mode='lines',
                                name=tip,line=dict(color='red',width=4),showlegend=False),)
        if np.isin('REAL', GEOL):
                Real.opacity=TGREAL
                fig.add_trace(Real)
                grealo={'display': 'block'}
        else:
            grealo={'display': 'none'}

        if np.isin('COL', GEOL):
                Colorado.opacity=TCOLORADO
                fig.add_trace(Colorado)
                coloradoo={'display': 'block'}
        else:
            coloradoo={'display': 'none'}       
        if np.isin('MUG', GEOL):
                Mugrosa.opacity=TMUGROSA
                fig.add_trace(Mugrosa)
                mugrosao={'display': 'block'}
        else:
            mugrosao={'display': 'none'}
        if np.isin('CHO', GEOL):
                Chorros.opacity=TCHORROS
                fig.add_trace(Chorros)
                chorroso={'display': 'block'}
        else:
            chorroso={'display': 'none'}
        if np.isin('EOC', GEOL):
                Eoceno.opacity=TEOCMED
                fig.add_trace(Eoceno)
                eocmedo={'display': 'block'}
        else:
            eocmedo={'display': 'none'}
        if np.isin('PER', CART):
                fig.add_trace(profile_plane(x0,y0,x1,y1))
        if np.isin('SEIS', PETRO):
                fig.add_trace(SISMICA)
                fig.update_layout(scene=dict(annotations=txts_p),
                overwrite=False)
        fig.update_layout(autosize=False,
                        width=850, height=600,
                        #margin=dict(l=50, r=50, b=50, t=50),
                        )
        fig.update_layout(
        scene = dict(aspectratio=dict(x=1,y=1.785714286,z=(42000/155540)*EXG),
                xaxis = dict(title='Longitud(°)',nticks=10, range=[loi,los]),
                yaxis = dict(title='Latitud(°)',nticks=10, range=[lai,las],),
                zaxis = dict(title='Elevación(msnm)',nticks=10, range=[-32000,10000],),),)

        return fig,START_DATE,grealo,coloradoo,mugrosao,chorroso,eocmedo,INYO

@app.callback(
     dash.dependencies.Output(component_id='Model_profile', component_property='figure'),



    [dash.dependencies.Input(component_id='DATE', component_property='start_date'),
     dash.dependencies.Input(component_id='DATE', component_property='end_date'),
     dash.dependencies.Input(component_id='MAGN', component_property='value'),
     dash.dependencies.Input(component_id='DEPTH', component_property='value'),
     dash.dependencies.Input(component_id='SEISMO', component_property='value'),
     dash.dependencies.Input(component_id='Longitud 1', component_property='value'),
     dash.dependencies.Input(component_id='Longitud 2', component_property='value'),
     dash.dependencies.Input(component_id='Latitud 1', component_property='value'),
     dash.dependencies.Input(component_id='Latitud 2', component_property='value') ])

def update_profile(START_DATE,END_DATE,MAGN,DEPTH,SEISMO,x0,x1,y0,y1):
    #Perfil sismico'----------------------------------------------------------------------------------------------------------------
    df_sismos_1=df_sismos[(df_sismos['FECHA - HORA UTC']<=END_DATE)&(df_sismos['FECHA - HORA UTC']>=START_DATE)&
    (df_sismos['MAGNITUD']>=MAGN[0])&(df_sismos['MAGNITUD']<=MAGN[1])
    &(df_sismos['PROF. (m)']>=DEPTH[0])&(df_sismos['PROF. (m)']<=DEPTH[1])]
    df_sismos_1['Unnamed: 0']=[i for i in range(0,len(df_sismos_1))]
    fig2= go.Figure()
    df_profile,dist_max=profile(x0,x1,y0,y1,df_sismos_1)
    text1=text_scatter(SEISMO,df_profile)
    profiles=go.Scatter(x=df_profile['DIST'],
                                y=df_profile['PROF. (m)'],
                                mode='markers',
                                name='Sismos',
                                error_y=dict(
                                    array=df_profile['ERROR PROFUNDIDAD SUP (m)'],                # set color to an array/list of desired values
                                    color='red',   # choose a colorscale
                                    symmetric=True,
                                    thickness=0.1,
                                    arrayminus=df_profile['ERROR PROFUNDIDAD (m)']
                                ),
                                marker=dict(size=df_profile['MAGNITUD']*5,
                                            color=df_profile['PROF. (m)'],
                                            colorscale='Jet',   # choose a colorscale
                                            opacity=0.8,
                                            cmax=100,
                                            cmin=-32000,),
                                            hovertemplate=text1)
    seismic_scale=go.Scatter(x=[dist_max/25,dist_max*1.2/25,dist_max*1.5/25,dist_max*2/25,dist_max*2.6/25],
                                y=[8000]*5,
                                 hoverinfo='text',
                                hovertext=np.array(['1','2','3','4','5']),
                                mode='markers',
                                name='Magnitudes',
                                showlegend=False,
                                marker=dict(size=(np.array([1,2,3,4,5])*5),
                                            opacity=0.8,
                                            ))
    for x_t,t_t in zip([dist_max/25,dist_max*1.2/25,dist_max*1.5/25,dist_max*2/25,dist_max*2.6/25],['1','2','3','4','5']):
        fig2.add_annotation(x=x_t, y=12200,hovertext=t_t,
            text=t_t,
            showarrow=False,
            yshift=0,
            
        font=dict(
            family="Courier New, monospace",
            size=12,
            
            )),
    fig2.add_trace(seismic_scale)  
    if dist_max>=2:
        scale_num=15;scale_text='15';scale_size=12
    elif dist_max<2 and dist_max>=1:
        scale_num=10;scale_text='10';scale_size=12
    else :
        scale_num=5;scale_text='5';scale_size=12
    fig2.add_annotation(x=(scale_num/111.1)/2, y=-38000,hovertext=scale_text+' km',
                text=scale_text+'km',
                showarrow=False,
                yshift=0,
            font=dict(
                family="Courier New, monospace",
                size=scale_size,
                
                ))
    fig2.add_trace(go.Scatter(x=np.linspace(0, scale_num/111.1, 4), y=np.array([-40000]*4),
                hovertext='Escala',
                mode='lines',
                name=scale_text+'km',
                showlegend=False,
                hoverinfo='text',
                marker=dict(color='black', size=8)))
    fig2.add_trace(go.Scatter(x=np.linspace(scale_num/111.1, scale_num*2/111.1, 4), y=np.array([-40000]*4),
                hovertext='Escala',
                mode='lines',
                hoverinfo='text',
                name=scale_text+'km',
                showlegend=False,
                marker=dict(color='white', size=8)))
    fig2.add_trace(go.Scatter(x=np.linspace(scale_num*2/111.1, scale_num*3/111.1, 4), y=np.array([-40000]*4),
                hovertext='Escala',
                mode='lines',
                hoverinfo='text',
                name=scale_text+'km',
                showlegend=False,
                marker=dict(color='black', size=8)))
    fig2.add_trace(go.Scatter(x=np.linspace(scale_num*3/111.1, scale_num*4/111.1, 4), y=np.array([-40000]*4),
                hovertext='Escala',
                hoverinfo='text',
                mode='lines',
                name=scale_text+'km',
                showlegend=False,
                marker=dict(color='white', size=8)))
    fig2.add_trace(go.Scatter(x=np.linspace(scale_num*4/111.1, scale_num*5/111.1, 4), y=np.array([-40000]*4),
                hovertext='Escala',
                hoverinfo='text',
                mode='lines',
                name=scale_text+'km',
                showlegend=False,
                marker=dict(color='black', size=8)))
    fig2.add_trace(profiles)      
    t_profile=topo_profile(x0,x1,y0,y1,df_topo)
    fig2.add_trace(t_profile)  
    Eoceno_1=geologic_profile(x0,y0,x1,y1,'datasets/DISCORDANCIA_EOCENO.txt','Discordancia del Eoceno Medio','red')
    Real_1=geologic_profile(x0,y0,x1,y1,'datasets/BASE_CUATERNARIO.txt','Tope Grupo Real','purple')
    Chorros_1=geologic_profile(x0,y0,x1,y1,'datasets/TOPE_CHORROS.txt','Tope Grupo Chorros','orangered')
    Colorado_1=geologic_profile(x0,y0,x1,y1,'datasets/TOPE_COLORADO.txt','Tope Formación Colorado','darkblue')
    Mugrosa_1=geologic_profile(x0,y0,x1,y1,'datasets/TOPE_MUGROSA.txt','Tope Formación Mugrosa','green')
    p1,p2=orientation(x0,y0,x1,y1)
    fig2.add_annotation(x=0, y=10000,
        text=p1,
        showarrow=False,
        yshift=0,
    font=dict(
        family="Courier New, monospace",
        size=25,
        ))
    fig2.add_annotation(x=dist_max, y=10000,
        text=p2,
        showarrow=False,
        font=dict(
        family="Courier New, monospace",
        size=25,
        ),
        yshift=0)
    fig2.add_traces(data=[Eoceno_1,Real_1,Chorros_1,Colorado_1,Mugrosa_1]) 
    fig2.update_layout(
                title="Perfil",
            ),
    fig2.update_layout(xaxis_title="Distancia (°)",
                        yaxis_title="Profundidad (m)")
    return fig2

@app.callback(
     dash.dependencies.Output(component_id='Iny_graph', component_property='figure'),



    [dash.dependencies.Input(component_id='TINY', component_property='value')])


def iny(TINY):
    datos_iny = pd.read_csv("datasets/inyeccion_geo.csv", delimiter = ';')
    name_campo=TINY
    fig = make_subplots(
        rows=1, cols=2,
        specs=[[{"type": "scatter"}, {"type": "bar"}]],
        subplot_titles=('Volumenes de agua (bbl) -'+name_campo, "Total por año" ))
    months=[]
    for i in datos_iny.columns:
        if '-' in i:
            months.append(i)
    name_campo='LA CIRA'
    iny_df=datos_iny[datos_iny['CAMPO']==name_campo]

    fig.add_trace(
        go.Scatter(x=months,y=[float(np.array(iny_df[x])[0]) for x in months],name=name_campo,showlegend=False),
        row=1, col=1
    )
    años=np.arange(2017,2022)
    old=np.array([0,0,0,0,0])
    for i in datos_iny['CAMPO']:
        if i!='TOTAL':
            campiny=datos_iny[datos_iny['CAMPO']==i]
            new=[float(np.array(campiny[age])[0]) for age in años.astype(str)]
            fig.add_trace(
                go.Bar(name=i, x=años, y=new+old,hovertemplate=['bbl:'+str(x) for x in new],showlegend=False),row=1, col=2)
            old=new
        else:
            pass
    fig.update_layout(barmode='stack')
    fig.update_yaxes(tickvals=np.arange(0,300000000+1,50000000),row=2, col=2)
    fig.update_xaxes(tickvals=años,row=2, col=2)
    fig.update_xaxes(tickangle=45)
    return fig

if __name__ == "__main__":
    app.run_server(debug=True)