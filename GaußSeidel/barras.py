# -*- coding: UTF-8 -*-
"""
Copyright (c) 2012 Itaman Cavalcanti de Oliveira <itamanman@yahoo.com.br>
Copyright (c) 2012 Tiago Luiz Satana de Souza    <tiagoluiz@oi.com.br>

Este programa é um software livre; você pode redistribui-lo e/ou 
modifica-lo dentro dos termos da Licença Pública Geral Menor GNU como 
publicada pela Fundação do Software Livre (FSF); na versão 3 da 
Licença, ou (na sua opnião) qualquer versão.

Este programa é distribuido na esperança que possa ser  util, 
mas SEM NENHUMA GARANTIA; sem uma garantia implicita de ADEQUAÇÂO
a qualquer MERCADO ou APLICAÇÃO EM PARTICULAR. Veja a Licença
Pública Geral GNU para maiores detalhes em <www.gnu.org>.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details at www.gnu.org
"""


import numpy
import math
import os
def tensao(adm_mat,barra,n_barras,q_reativa,p_ativa,v):
    soma=0.0
    for cnt in range(n_barras):
        if cnt!=barra:
            soma=soma+adm_mat[barra,cnt]*v[cnt]
    v_gauss=(1/adm_mat[barra,barra])*(((p_ativa[barra]-(1j)*(q_reativa[barra]))/(v[barra].conjugate()))-soma)
    return v_gauss

def reativa(adm_mat,barra,n_barras,v):
    soma2=0.0
    for cont in range(n_barras):
        soma2=soma2+(v[barra].conjugate())*((v[cont])*adm_mat[barra,cont])
    q_gauss=(-1)*soma2.imag
    return q_gauss

def ativa(adm_mat,barra,n_barras,v):
    soma2=0.0
    for cont in range(n_barras):
        soma2=soma2+(v[barra].conjugate())*((v[cont])*adm_mat[barra,cont])
    p_gauss=soma2.real
    return p_gauss


def mudar_ele_mat():
    bur=raw_input('\nHouve algum erro na entrada de elementos? Digite (s ou n): ')
    if bur=='s':
        print 'Mude o elemento da matriz que há o erro'
        os.system("kwrite vetor_impe_li.txt")
        bu=numpy.loadtxt("vetor_impe_li.txt", dtype="|S25", delimiter=",", skiprows=0)
        print bu
        cont3=0
        adm_mat2=numpy.array([[0]*n_barras]*n_barras, complex)
        
        for i in range(n_barras):
            for k in range(n_barras):
                if i<k:
                    if complex(bu[cont3])==0:
                        adm_mat2[i,k]=0
                        adm_mat2[k,i]=adm_mat2[i,k]
                    else:
                        adm_mat2[i,k]=1/complex(bu[cont3])
                        adm_mat2[k,i]=adm_mat2[i,k]
                    cont3=cont3+1
        for j in range(n_barras):
            for k in range(n_barras):
                if j!=k:
                    adm_mat2[j,j]= adm_mat2[j,j]+adm_mat2[j,k]
                    adm_mat[j,k]=(-1)*adm_mat2[j,k]
            adm_mat[j,j]= adm_mat2[j,j]+shunt[j]
        print '\nMatriz Admintancia:\n', adm_mat




'''def inversa_mat_adm(adm_mat):
    z_bus=numpy.linalg.inv(adm_mat)
    return z_bus'''
   
def falta(barra_f,adm_mat,v):
    z_f=complex(input('\nImpedância de falta: '))
    z_bus=numpy.linalg.inv(adm_mat)
    iff=v[barra_f]/(z_f+z_bus[barra_f,barra_f])
    for i in range(n_barras):
        v_f[i]=v[i]-(z_bus[i,barra_f]/(z_bus[barra_f,barra_f]+z_f))*v[barra_f]
    for j in range(n_barras):
        for k in range(n_barras):
            if j!=k:
                i_l[j,k]=adm_mat[j,k]*(v_f[j]-v_f[k])
    print 'Corrente de falta: ', iff
    for j in range(n_barras):
        for k in range(n_barras):
            jj=j+1
            kk=k+1
            if j<k:
                print 'Corrente de falta entre as barras %d-%d: ' %(jj,kk), i_l[j,k]
    for i in range(n_barras):
        ii=i+1
        print 'Tensão pós falta na barra %d: ' %ii, v_f[i]
'''
oi
    z_f=complex(input('\nImpedância de falta: '))
    print adm_mat,'\n'
    z_bus=numpy.linalg.inv(adm_mat)
    print z_bus
    If=v[barra_f]/(z_f+z_bus[barra_f,barra_f])
    for i in range(n_barras):
        v_f[i]=v[i]-(z_bus[i,barra_f]/(z_bus[barra_f,barra_f]+z_f))*v[barra_f]
    for j in range(n_barras):
        for k in range(n_barras):
            if j!=k:
                i_l[j,k]=adm_mat[j,k]*(v_f[j]-v_f[k])
    print 'Corrente de falta: ', If
    for j in range(n_barras):
        for k in range(n_barras):
            jj=i+1
            kk=k+1
            print 'Corrente de falta entre as barras %d-%d' %(jj,kk), i_l[j,k]
    for i in range(n_barras):
        print 'Tensão pos falta: ', v_f[i]'''
    
    
def carregar_n_mat():
    print '\n'
    for i in range(n_barras):
        ii=i+1
        sh=raw_input('Impedância shunt da Barra %d: ' %ii)
        if sh=='ne' or sh=='0':
            shunt[i]=0
        else:    
            shunt[i]=1/complex(sh)
    print '\n'
    cont2=0
    for j in range(n_barras):
        for k in range(n_barras):
            t=k+1
            u=j+1
            if j<k:
                impedancia=raw_input('Impedância da linha %d-%d: ' %(u,t))
                if impedancia=='ne':
                    z_impe[cont2]=0
                    adm_mat[j,k]=0
                    adm_mat[k,j]=adm_mat[j,k]
                else:
                    z_impe[cont2]=complex(impedancia)
                    adm_mat[j,k]=1/complex(impedancia)
                    adm_mat[k,j]=adm_mat[j,k]
                cont2=cont2+1
    for j in range(n_barras):
        for k in range(n_barras):
            if j!=k:
                adm_mat[j,j]= adm_mat[j,j]+adm_mat[j,k]
                adm_mat[j,k]=(-1)*adm_mat[j,k]
        adm_mat[j,j]=adm_mat[j,j]+shunt[j]
    print '\nMatriz Admintância:\n', adm_mat
    numpy.savetxt("vetor_impe_li.txt", z_impe, fmt="%10.8f", delimiter=",")
    mudar_ele_mat()

def tipo_barras(h):
    print '\nConsidere a Barra de Referência como a Barra 1'
    v_mod=input('\nMódulo da Tensão da Barra de Referência (em p.u.): ')
    v_fase=input('Fase da Tensão da Barra de Referência (em Graus): ')
    v[0]=v_mod*(math.cos((v_fase*math.pi)/180)+1j*math.sin((v_fase*math.pi)/180))

    for i in range(n_barras):
        ii=i+1
        if i==0:
            tipo[i]==['slack']
        else:
            tipo[i]=[raw_input('\nTipo da barra %d: ' %ii)]
            while tipo[i] not in[['pq'], ['pv'],['PV'],['PQ']]:
                tipo[i]=[raw_input('Tipo da barra %d: ' %ii)]
        if tipo[i]==['pq'] or tipo[i]==['PQ']:
            p_ativa[i]=input('\nPotência Ativa da Barra %d (em p.u.): ' %ii)
            q_reativa[i]=input('Potência Reativa da Barra %d:(em p.u.): ' %ii)
        if tipo[i]==['pv'] or tipo[i]==['PV']:
            p_ativa[i]=input('\nPotência Ativa da Barra %d (em p.u.): ' %ii)
            v_mod=input('Módulo da Tensão da Barra %d (em p.u.): ' %ii)
            v_fase=0 # condicao inicial
            v[i]=v_mod*(math.cos((v_fase*math.pi)/180)+1j*math.sin((v_fase*math.pi)/180))
    numpy.savetxt("vetor_tensao.txt", v, fmt="%10.8f", delimiter=",")
    numpy.savetxt("vetor_p.txt", p_ativa, fmt="%10.8f", delimiter=",")
    numpy.savetxt("vetor_q.txt", q_reativa, fmt="%10.8f", delimiter=",")
    if h==1:
        carregar_n_mat()    

def mudar_tensao():
    bur2=raw_input('\nHouve algum erro nas tensões? (s ou n): ')
    while bur2 not in ['s', 'n']:
        bur2=raw_input('(s ou n): ')
    if bur2=='s':
        print "Mude os valores:"
        os.system("kwrite vetor_tensao.txt")
        bu2=numpy.loadtxt("vetor_tensao.txt", dtype="|S25", delimiter=",",skiprows=0)
        for i in range(n_barras):
            v[i]=complex(bu2[i])
        
def mudar_potencias():
    bur3=raw_input('\nHouve algum erro nas potências? (s ou n): ')
    while bur3 not in ['s', 'n']:
        bur3=raw_input('(s ou n): ')
    if bur3=='s':
        print "Mude os valores das potências ativas:"
        os.system("kwrite vetor_p.txt")
        bu3=numpy.loadtxt("vetor_p.txt", dtype="|S25", delimiter=",",skiprows=0)
        for i in range(n_barras):
            p_ativa[i]=float(bu3[i])
        print "Mude os valores das potências reativas:"
        os.system("kwrite vetor_q.txt")
        bu4=numpy.loadtxt("vetor_p.txt", dtype="|S25", delimiter=",",skiprows=0)
        for i in range(n_barras):
            q_reativa[i]=float(bu4[i])
        print p_ativa, ' ', q_reativa
    
def gauss_seidel():
    erro=0.00001 #Erro pr
    erro1=0.0
    erro2=0.0
    erro3=1.0
    if n_barras>1:
        cont=0
        while (erro3 > erro):
            cont=cont+1
            for i in range(n_barras):
                if tipo[i]==['pv']:
                    barra=i
                    q_reativa[i]=reativa(adm_mat,barra,n_barras,v)
                    v_cal=tensao(adm_mat,barra,n_barras,q_reativa,p_ativa,v)
                    tan_v=v_cal.imag/v_cal.real
                    ang_c=math.atan(tan_v)
                    v_ant=v[i]
                    '''v[i]=abs(v_ant)*(math.cos(math.atan(ang_c))+1j*math.sin(math.atan(ang_c)))'''
                    v[i]=abs(v_ant)*(math.cos(ang_c)+1j*math.sin(ang_c))
                    erro1=abs(v[i]-v_ant)

                if tipo[i]==['pq']:
                    barra=i
                    v_ant2=v[i]
                    v[i]=tensao(adm_mat,barra,n_barras,q_reativa,p_ativa,v)
                    erro2=abs(v[i]-v_ant2)
            if erro1>erro2:
                erro3=erro1
            if erro2>erro1:
                erro3=erro2
    for i in range(n_barras):
        ii=i+1
        tan=v[i].imag/v[i].real
        ang=(180*math.atan(tan))/math.pi
        print '\nTensao na barra %d :'%ii, v[i]
        print 'Módulo: ', abs(v[i]), 'Ângulo de fase (em graus): ', ang
    
    print '\nPotência Barra de Referência', 'P: ', ativa(adm_mat,0,n_barras,v), 'Q: ',reativa(adm_mat,0,n_barras,v)
    print '\nQuantidade de Iterações:', cont
    
open('vetor_p','w')
open('vetor_q','w')
open('vetor_tensao.txt','w')
open('vetor_impe_li.txt','w')
inicio=''
while inicio=='':
    inicio=raw_input('\nIniciar programa (digite I)\nSair (digite S): ')
    while inicio not in ['I','s', 'S','i']:
        inicio=raw_input('Iniciar programa (digite I)\nSair (digite S): ')
    if inicio in ['i','I']:
        inicio2=raw_input('\nNovo sitema? (s ou n): ')
        while inicio2 not in ['s', 'n', 'S', 'N']:
            inicio2=raw_input('Novo sitema? (s ou n) :')
        if inicio2 in ['s', 'S']:
            n_barras = int(input('\nDigite a quantidade de barras do sistema: '))
            tipo=['']*n_barras
            j=0
            for i in range(n_barras):
                j=i+j
            z_impe=numpy.array([0]*j, complex)
            impe=[[0]*n_barras]*n_barras
            z_bus=numpy.array([0]*j, complex)
            adm_mat=numpy.array(impe, complex)
            adm_mat2=numpy.array(impe, complex)
            p_ativa=numpy.array([0]*n_barras, float)
            q_reativa=numpy.array([0]*n_barras, float)
            v=numpy.array([1+0j]*n_barras, complex)
            v_f=numpy.array([1+0j]*n_barras, complex)         #condicao inicial V's=0
            i_l=numpy.array([[0]*n_barras]*n_barras, complex)
            shunt=[0]*n_barras
            tipo_barras(1)
            gauss_seidel()
            raw_input('\nAperte ENTER para voltar a tela inicial')
            inicio=''
            os.system("clear")
        if inicio2 in ['n','N']:
            passo3=raw_input('\nMudar matriz admitância (M)\nMudar tipos de Barras (B)\nMudar Tensões (T)\nMudar Potências (P)\n Calcular Falta (F): ')
            while passo3 not in ['m', 'b', 'B', 'M', 'P','p','T','t', 'f','F']:
                passo3=raw_input('Mudar matriz admitancia (M)\nMudar tipos de Barras (B): ')
                
            if passo3 in ['f','F']:
                barra_ch=input('Barra de falta: ')
                barra_f=barra_ch-1
                falta(barra_f,adm_mat,v)
            if passo3 in ['m', 'M']:
                passo4=raw_input('\nCarregar nova Matriz (C)\nMudar elemento (M): ')
                while passo4 not in ['m', 'c', 'C', 'M']:
                    passo4=raw_input('Carregar nova Matriz (C)\nMudar elemento (M): ')
                if passo4 in ['c', 'C']:
                    adm_mat=numpy.array(impe, complex)
                    carregar_n_mat()
                if passo4 in ['m', 'M']:
                    mudar_ele_mat()
            if passo3 in ['b', 'B']:
                tipo_barras(0)
            if passo3 in ['t', 'T']:
                mudar_tensao()
            if passo3 in ['p','P']:
                mudar_potencias()
            gauss_seidel()
            raw_input('\nAperte ENTER para voltar a tela inicial')
            inicio=''
            os.system("clear")
    if inicio in ['s','S']:
        open('vetor_impe_li.txt').close()
        open('vetor_p').close()
        open('vetor_q').close()
        open('vetor_tensao.txt','w').close()
        print 'Fechando...'
        exit