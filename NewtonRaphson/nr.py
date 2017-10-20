# -*- coding: UTF-8 -*-
"""
Copyright (c) 2013 Itaman Cavalcanti de Oliveira <itamanman@yahoo.com.br>

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
from libsvn import ra

'''def fluxo_p(n_barras,v,adm_mat,f_p):
    p=f_p
    for i in range(n_barras):
        for j in range(n_barras):
            p[i,j]=-abs(v[i])*abs(adm_mat[i,j])*math.cos(numpy.angle(adm_mat[i,j]))+abs(v[i])*abs(v[j])*abs(adm_mat[i,j])*math.cos(-numpy.angle(adm_mat[i,j])+numpy.angle(v[i])-numpy.angle(v[j]))
    return p

def fluxo_q(n_barras,v,adm_mat,f_q):
    q=f_q
    for i in range(n_barras):
        for j in range(n_barras):
            q[i,j]=-abs(v[i])*abs(adm_mat[i,j])*math.sin(numpy.angle(adm_mat[i,j]))+abs(v[i])*abs(v[j])*abs(adm_mat[i,j])*math.sin(-numpy.angle(adm_mat[i,j])+numpy.angle(v[i])-numpy.angle(v[j]))
    return q

def fl(v, adm_mat):
    f=0
    for i in range(nb):
        for j in range(nb):
            if i!=j:
                f=v[i]*(v[i]-v[j])/adm_mat[i,j]'''

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
        print '\nMatriz Admintância:\n', adm_mat


def p_NR(barra, n_barras, p_ativa, v, adm_mat):
    soma_p=0.0
    for i in range(n_barras):
        soma_p=soma_p+abs(v[i])*abs(adm_mat[barra,i])*math.cos(-numpy.angle(adm_mat[barra,i])+numpy.angle(v[barra])-numpy.angle(v[i]))
    p_n=p_ativa[barra]-abs(v[barra])*soma_p
    return p_n
def q_NR(barra, n_barras, q_reativa, v, adm_mat):
    soma_q=0.0
    for i in range(n_barras):
        soma_q=soma_q+abs(v[i])*abs(adm_mat[barra,i])*math.sin(-numpy.angle(adm_mat[barra,i])+numpy.angle(v[barra])-numpy.angle(v[i]))
    
    q=q_reativa[barra]-abs(v[barra])*soma_q
    return q
def jacobiano(n_barras, adm_mat, v,n,jaco,tipo,n_pv):
    nb=n_barras
    jacob=jaco
    if n_pv==0:
	for i in range(1,nb):
	    for j in range(1,nb):
		if i!=j:
		    jacob[i-1,j-1]=abs(v[i])*abs(v[j])*abs(adm_mat[i,j])*math.sin(numpy.angle(v[i])-numpy.angle(v[j])-numpy.angle(adm_mat[i,j]))
		if i==j:
		    soma_j1=0.0
		    for j1 in range(n_barras):
			if i!=j1:
			    soma_j1=soma_j1+abs(adm_mat[i,j1])*abs(v[j1])*math.sin(numpy.angle(v[i])-numpy.angle(v[j1])-numpy.angle(adm_mat[i,j1]))
		    jacob[i-1,j-1]=-abs(v[i])*soma_j1
	for i in range(1,nb):
	    for j in range(1,nb-n_pv):
		if i!=j:
		    jacob[i-1,j-nb+n]=abs(v[i])*abs(adm_mat[i,j])*math.cos(numpy.angle(v[i])-numpy.angle(v[j])-numpy.angle(adm_mat[i,j]))
		if i==j:
		    soma_j2=0.0
		    for j2 in range(n_barras):
			soma_j2=soma_j2+abs(v[j2])*abs(adm_mat[i,j2])*math.cos(numpy.angle(v[i])-numpy.angle(v[j2])-numpy.angle(adm_mat[i,j2]))
		    jacob[i-1,j-nb+n]=soma_j2+abs(v[i])*abs(adm_mat[i,i])*math.cos(numpy.angle(adm_mat[i,i]))
	for i in range(1,nb-n_pv):
	    for j in range(1,nb):
		if i!=j:
		    jacob[i-nb+n,j-1]=-abs(v[i])*abs(v[j])*abs(adm_mat[i,j])*math.cos(numpy.angle(v[i])-numpy.angle(v[j])-numpy.angle(adm_mat[i,j]))
		if i==j:
		    soma_j3=0.0
		    for j3 in range(n_barras):
			if i!=j3:
			    soma_j3=soma_j3+abs(adm_mat[i,j3])*abs(v[j3])*math.cos(numpy.angle(v[i])-numpy.angle(v[j3])-numpy.angle(adm_mat[i,j3]))
		    jacob[i-nb+n,j-1]=abs(v[i])*soma_j3
	for i in range(1,nb-n_pv):
	    for j in range(1,nb-n_pv):
		if i!=j:
		    jacob[i-nb+n,j-nb+n]=abs(v[i])*abs(adm_mat[i,j])*math.sin(numpy.angle(v[i])-numpy.angle(v[j])-numpy.angle(adm_mat[i,j]))
		if i==j:
		    soma_j4=0.0
		    for j4 in range(nb):
			soma_j4=soma_j4+abs(v[j4])*abs(adm_mat[i,j4])*math.sin(numpy.angle(v[i])-numpy.angle(v[j4])-numpy.angle(adm_mat[i,j4]))
		    jacob[i-nb+n,j-nb+n]=soma_j4-abs(v[i])*abs(adm_mat[i,i])*math.sin(numpy.angle(adm_mat[i,i]))
    else:
	for i in range(1,nb):
	    for j in range(1,nb):
		if i!=j:
		    jacob[i-1,j-1]=abs(v[i])*abs(v[j])*abs(adm_mat[i,j])*math.sin(numpy.angle(v[i])-numpy.angle(v[j])-numpy.angle(adm_mat[i,j]))
		if i==j:
		    soma_j1=0.0
		    for j1 in range(n_barras):
			if i!=j1:
			    soma_j1=soma_j1+abs(adm_mat[i,j1])*abs(v[j1])*math.sin(numpy.angle(v[i])-numpy.angle(v[j1])-numpy.angle(adm_mat[i,j1]))
		    jacob[i-1,j-1]=-abs(v[i])*soma_j1
	for i in range(1,nb):
	    for j in range(1,nb-n_pv):
		if i!=j:
		    jacob[i-1,j-(n_pq+n_pv)+n]=abs(v[i])*abs(adm_mat[i,j])*math.cos(numpy.angle(v[i])-numpy.angle(v[j])-numpy.angle(adm_mat[i,j]))
		if i==j:
		    soma_j2=0.0
		    for j2 in range(n_barras):
			soma_j2=soma_j2+abs(v[j2])*abs(adm_mat[i,j2])*math.cos(numpy.angle(v[i])-numpy.angle(v[j2])-numpy.angle(adm_mat[i,j2]))
		    jacob[i-1,j-(n_pq+n_pv)+n]=soma_j2+abs(v[i])*abs(adm_mat[i,i])*math.cos(numpy.angle(adm_mat[i,i]))
	for i in range(1,nb-n_pv):
	    for j in range(1,nb):
		if i!=j:
		    jacob[i-(n_pq+n_pv)+n,j-1]=-abs(v[i])*abs(v[j])*abs(adm_mat[i,j])*math.cos(numpy.angle(v[i])-numpy.angle(v[j])-numpy.angle(adm_mat[i,j]))
		if i==j:
		    soma_j3=0.0
		    for j3 in range(n_barras):
			if i!=j3:
			    soma_j3=soma_j3+abs(adm_mat[i,j3])*abs(v[j3])*math.cos(numpy.angle(v[i])-numpy.angle(v[j3])-numpy.angle(adm_mat[i,j3]))
		    jacob[i-(n_pq+n_pv)+n,j-1]=abs(v[i])*soma_j3
	for i in range(1,nb-n_pv):
	    for j in range(1,nb-n_pv):
		if i!=j:
		    jacob[i-(n_pq+n_pv)+n,j-(n_pq+n_pv)+n]=abs(v[i])*abs(adm_mat[i,j])*math.sin(numpy.angle(v[i])-numpy.angle(v[j])-numpy.angle(adm_mat[i,j]))
		if i==j:
		    soma_j4=0.0
		    for j4 in range(nb):
			soma_j4=soma_j4+abs(v[j4])*abs(adm_mat[i,j4])*math.sin(numpy.angle(v[i])-numpy.angle(v[j4])-numpy.angle(adm_mat[i,j4]))
		    jacob[i-(n_pq+n_pv)+n,j-(n_pq+n_pv)+n]=soma_j4-abs(v[i])*abs(adm_mat[i,i])*math.sin(numpy.angle(adm_mat[i,i]))
    return jacob

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
    
def newton_raphson():
    erro=0.00001 #Erro pr
    erro3=1.0
    if n_barras>1:
        cont=0
        while (erro3 > erro):
            for i in range(n):
                for j in range(n):
                    if abs(delta_p[i])>abs(delta_p[j]):
                        erro3=abs(delta_p[i])
                    if abs(delta_p[j])>abs(delta_p[i]):
                        erro3=abs(delta_p[j])
            cont=cont+1
            contador=n_pq+n_pv
            for i in range(n_barras):
                barra=i
                if tipo[i] in [['PQ'],['pq']]:
                    delta_p[i-1]=p_NR(barra, n_barras, p_ativa, v, adm_mat)
                    delta_p[contador]=q_NR(barra, n_barras, q_reativa, v, adm_mat)
                    contador=contador+1     
                if tipo[i] in [['PV'],['pv']]:
                    delta_p[i-1]=p_NR(barra, n_barras, p_ativa, v, adm_mat)
        
            
            jacobia=jacobiano(n_barras, adm_mat, v, n, jaco,tipo, n_pv)
            delta_re=numpy.dot(numpy.linalg.inv(jacobia),delta_p)
            
            for i in range(n_pv+n_pq):
                delta_ang[i]=delta_re[i]
            
            for i in range(n_pq+n_pv,((2*n_pq)+n_pv)):
                delta_ten[i-(n_pq+n_pv)]=delta_re[i]
            
            counter=0
            for i in range(n_barras):
                
                if tipo[i] in [['PQ'], ['pq']]:
                    v_ant=v[i]
                    v[i]=(abs(v_ant)+delta_ten[counter])*(math.cos(numpy.angle(v_ant)+delta_ang[i-1])+1j*math.sin(numpy.angle(v_ant)+delta_ang[i-1]))
                    counter=counter+1
                if tipo[i] in [['pv'], ['PV']]:
                    v_ant=v[i]
                    v[i]=abs(v_ant)*(math.cos(numpy.angle(v_ant)+delta_ang[i-1])+1j*math.sin(numpy.angle(v_ant)+delta_ang[i-1]))
                
                       
    for i in range(n_barras):
        ii=i+1
        tan=v[i].imag/v[i].real
        ang=(180*math.atan(tan))/math.pi
        print '\nTensao na barra %d :'%ii, v[i]
        print 'Módulo: ', abs(v[i]), 'Ângulo de fase (em graus): ', ang
    
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
            flu=numpy.array(impe, complex)
            f_p=numpy.array(impe, complex)
            f_q=numpy.array(impe, complex)
            adm_mat2=numpy.array(impe, complex)
            p_ativa=numpy.array([0]*n_barras, float)
            q_reativa=numpy.array([0]*n_barras, float)
            v=numpy.array([1+0j]*n_barras, complex)
            v_f=numpy.array([1+0j]*n_barras, complex)         #condicao inicial V's=0
            i_l=numpy.array([[0]*n_barras]*n_barras, complex)
            shunt=[0]*n_barras
            tipo_barras(1)
            n_pq=0
            n_pv=0
            for l in range(n_barras):
                if tipo[l] == ["pq"]:
                    n_pq = n_pq + 1
                if tipo[l] == ["pv"]:
                    n_pv = n_pv + 1
            n=n_pv+2*n_pq
            jaco=numpy.array([[0.0]*n]*n)
            delta_p=numpy.array([1.0]*n)
            delta_ang=numpy.array([0.0]*(n_pq+n_pv))
            delta_ten=numpy.array([0.0]*n_pq)   
            newton_raphson()
            raw_input('\nAperte ENTER para voltar a tela inicial')
            inicio=''
            os.system("clear")
        if inicio2 in ['n','N']:
            passo3=raw_input('\nMudar matriz admitância (M)\nMudar tipos de Barras (B)\nMudar Tensões (T)\nMudar Potências (P): ')
            while passo3 not in ['m', 'b', 'B', 'M', 'P','p','T','t', 'f','F']:
                passo3=raw_input('Mudar matriz admitancia (M)\nMudar tipos de Barras (B): ')
                
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
            newton_raphson()
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