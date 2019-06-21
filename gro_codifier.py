import string, pprint, random, itertools
import numpy as np


def designer(n_species,n_interactions,av_force, n_matrixes, n_pos):
#El input es el n de especies que va a simular, el n de interacciones y la fuerza media, que sera la desviación típica de la distribución
#Distribución normal de la cual extraer las fuerzas de cada interacción. También el número de matrices, que sera el n de archivos .gro

    n_zeros=n_species*(n_species-1)-n_interactions                  #Matriz n*n-1 porque no queremos la diagonal, realmente es un vector
    zeros = [0 for x in range(n_zeros)]                             #Esa matriz es realmente un vector zeros + non-zeros, estos ultimos salen de una normal

    ##no zeros ahora son el mismo valor con distintos signos
    ##n_pos es el n de interacciones con signo positivo
    nonz_neg = [av_force for x in range(n_interactions-n_pos)]
    nonz_pos = [-av_force for x in range(n_pos)]
    prematrix = zeros+nonz_neg+nonz_pos
    prematrixes_list=[]                                             #Vamos a generar algunas matrices posibles, para ello mezclamos los elementos del la lista matriz
    for i in range(n_matrixes):
        random.shuffle(prematrix)                                   #Generamos algunas posibles prematrices por aleatorización. Podría salir repetidas pero es muy improbable. Se podrían hacer permutaciones pero es inviable computacionalmente (Memory Error).
        prematrixes_list.append(list(prematrix))

    matrixes_list=[]
    for element in prematrixes_list:                                #Modificamos cada vector/lista para añadir los 0 de la diagonal

            for i in range(0,n_species**2,n_species+1):             #Añadir ceros a la diagonal de la matriz (realmente una lista/vector) 
                    element.insert(i,0)
            #pprint.pprint(element)#QUITAR
            matrixes_list.append(element)                           #Se obtiene el conjunto de las matrices completas. tamaño=n_matrices
            
    output_dicts=[]
    for matrix in matrixes_list:
            species_dict={}
            letters=list(string.ascii_uppercase)                    #Esto permite nombrar luego cada metabolito de interacción con una letra mayúscula
            for i in range(1,n_species+1):
                    species_dict[i]={"out":[],"in":[],"force":[]}   #Cada especie tiene una serie de metabolitos que produce(out) y que absorbe(in)

            for j in range(1,n_species+1):                          #Con cada matriz se procede de este modo para elaborar el esquema metabolico de cada especie
                    for k in range(1,n_species+1):                  #Las variables j y k son los indices de la matriz con los que establecer el metabolismo de cada especie
                            value=matrix.pop(0)                     #Extraemos el valor correpondiente y si no es 0 se establece interacción y se modifica el metabolismo de las especies corresponientes
                            if value!=0:
                                    metabolite=letters.pop(0)                   #Definimos un metabolito con el instanciador para que medie la interaccion
                                    species_dict[j]["out"].append(metabolite)   #Despues asignamos este metabolito a la lista correspondienta para cada especie
                                    species_dict[k]["in"].append(metabolite)
                                    species_dict[k]["force"].append(value)#Se asigna el valor en la fuerza de interacción que modificara la tasa de crecimiento

            output_dicts.append(species_dict)
    pprint.pprint(output_dicts)
    return output_dicts                                             #Se devuelve el conjunto de diccionarios (esquemas metabólicos), con los que hacer los archivos .gro 

#ejemplo=encoder(3,3,1,3)

#pprint.pprint(ejemplo)

midline_n = 56                                                      #Asignamos una linea de división para el modelo que vamos a leer
model_file = open ("modelo2.gro", "r")                               #Leemos el modelo del cual sacar el comienzo y el final del archivo .gro
model = model_file.readlines()
model_file.close()

model_init = model[:midline_n]                                      #Usando la linea umbral separamos entre lines de comienzo y lineas del final de un .gro
model_end = model[midline_n:]

def encoder(file, scheme):

    for species in scheme.keys():                                   #Cada especie esta asociada a un diccionario, esquema de su metabolismo

        metabolism = scheme[species]

        #Comenzando por los metabolitos en out, por cada uno hay que codificar una señal
        for signal in metabolism ['out']:
            ####SEÑALES
            molecule=signal+' := signal([ name := "'+signal+'", kdiff := diff_k, kdeg := degrad_k ]);\n\n'
            file.write(molecule)
            
        #Ahora se escribe lapuerta (metabolica) el operon y plásmido para cada especie. Y LA PERTURBACIÓN!!!!!!!!!!!!
        gate='molecule_gate([ name := "gate'+str(species)+'",\n\tfunction := "YES",\n\tinput := {"m'+str(species)+'"}\n]);\n\n'
        file.write(gate)

    file.write('\n//Parte metabolica\n\n')

    for species in scheme.keys():
        name_line='metabolism ([name := "met'+str(species)+'",\n'
        gate_line='\t\t\tgate := "gate'+str(species)+'",\n'                     #IMPORTANTE
        metabolites=''
        dev_line='\t\t\tdeviations := {'
        deviations='0.0},\n'
        morelines='\t\t\tnoise_time := 20.0,\n\t\t\tsensitivity_steps := 2.0,\n\n'

        metabolism = scheme[species]
        fluxes='"biomass"'
        fluxes_in=''
        fluxes_w=''
        fluxes_lines=[]
        for metabolite in metabolism['out']:
            metabolites='"'+metabolite+'", '+metabolites
            flux=metabolite+'out'
            fluxes='"'+flux+'", '+fluxes
            flux_lines='\t\t\t'+flux+' := [metabolite := "'+metabolite+'",\n\t\t\t\t\t f := [ bias := 0.1 ]\n\t\t\t],\n'
            fluxes_lines.append(flux_lines)

        for metabolite in metabolism['in']:
            metabolites='"'+metabolite+'", '+metabolites            #Rellenar flujos
            flux=metabolite+'in'                                    #Nombrar flujos
            fluxes_in='"'+flux+'", '+fluxes_in 
            fluxes='"'+flux+'", '+fluxes                            #Rellenar metabolitos de entrada. Hay un umbral para que comience a absorber.
            flux_lines='\t\t\t'+flux+' := [metabolite := "'+metabolite+'",\n\t\t\t\t\t functions := {"ins", "max"},\n\t\t\t\t\t ins := [ metabolites := {"'+metabolite+'"}, metabolites_w := {-0.5}, bias := 0.0 ],\n\t\t\t\t\t max := [ bias := -0.02 ],\n\t\t\t\t\t tree := [ '+metabolite+' := threshold_up, up := "max", low := "ins"]\n\t\t\t],\n'
            fluxes_lines.append(flux_lines)

       # if fluxes_lines:
        #    fluxes_lines[-1]=fluxes_lines[-1][:-2]

        for w in metabolism['force']:                               #No es necesario el redondeo de las fuerzas de interacción
            w=round(w,3)
            if w>0:
                fluxes_w='f_pos'+', '+fluxes_w
            else:
                fluxes_w='f_neg'+', '+fluxes_w

        if metabolites:
            metabolites=metabolites[:-2]
        if fluxes_in:
            fluxes_in=fluxes_in[:-2]
            fluxes_w=fluxes_w[:-2]
            biomass_line='\t\t\tbiomass := [ metabolite := "biomass",\n\t\t\t\t\t\tf := [ fluxes := {'+fluxes_in+'}, fluxes_w := {'+fluxes_w+'}, bias := growth_rate ]\n\t\t\t]\n]);\n\n'

        else:
            biomass_line='\t\t\tbiomass := [ metabolite := "biomass",\n\t\t\t\t\t\tf := [ bias := growth_rate ]\n\t\t\t]\n]);\n\n'


        metab_line='\t\t\tmetabolites := {'+metabolites+'},\n'
        flux_line='\t\t\tfluxes := {'+fluxes+'},\n'
        dev_line=dev_line+deviations

        #biomass_line=']);\n\n'
        ###YA NO SE PONE NI DEV NI EL RESTO, METAB SOLO SI HAY UN TREE, POR AHORA NO FUNCIONA
        file.writelines([name_line,gate_line, metab_line,flux_line])
        file.writelines(fluxes_lines+[biomass_line])
        
    return


def programer(species_schemes):                                     #Se introduce el conjunto de diccionarios (esquemas)

    for i in range(len(species_schemes)):                           #Cada esquema es un diccionario con n especies

        scheme = species_schemes[i]
        experiment_n="exp"+str(i)
        file = open(experiment_n+".gro", "w")                       #Se genera un archivo .gro para cada diseño

        model_endmod=list(model_end)                                #Generamos una nueva lista (evitar el ) para añadir el fopen y escribirlo después
        outputfile='fp := fopen (route1 <> "'+experiment_n+'.dat", "w" );\n'
        model_endmod.insert(15,outputfile)
        #Escritura de archivos
        file.writelines(model_init)
        encoder(file, scheme)                                       #Y se escribe en el codigo correspondiente al diseño introducido
        file.writelines (model_endmod)
        file.close()

    return

    

    
programer(designer(5,14,0.5,10,14))
                        











