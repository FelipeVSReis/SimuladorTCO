                                                                        /*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package SimuladorTCO;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;                                                            

/**
 *
 * @author Casa
 */
public class SimuladorTCO_fotoMod{
    //parametros do modelo manhatttan
    int N=14; //Número de Blocos numa fileira
    int n=12; //Número de Resisdencias numa fileira
    double l=0.03; // Dsitancia entre residencias (km)
    double L=n*l;        // Distancia Entre Blocos
    double num_andares=1;
    //parametros de custos
    double sal=95; //($/h)
    double PrKw=0.58; // Preço Energia
    double prAluguel= 183.96;//Preco aluguel por ano
    double prGerenc = 12.300*4*12; // Preco de gerenciamento que acabou nao sendo usado rs
    double cAluguel = 0; // Custo com aluguel
    double cGerenc = 0; // Custo com Gerenciamento
    //double FatorMultiplicador=10;
    double vel=20; //(Km/h)
    int SR=8;   //splitting ratio(16,32 ou 64)
    int NPC=72;  //numero de portas olt por chassi botar 72
   //subtituir e arqpad esqprot
   //seletores
    int typePon = 0;
    // 0 GPON
    // 1 10GPON
    // 2 40GPON
    int typeSF = 1;
    // 0 sem painel
    // 1 com painel
    //Parametros Monte Carlo
    //int numittt=200;
    int Ny=5;
    int Num_Tentativas=1000;
    double PorcRep=0.3; //reparo
    //Dados de Equipamentos
    // {onu,fiber,splitter,olt port,Rn chassi,Olt Chassi,switches,GES,Micro,Pico,Antena,macro}
    //PICO, ONU, Painel ONU, Inversor ONU, Distribution Fiber Step, SPL, RN Chassi, Feeder Fiber Step, OLT, OltChassi, Painel OLT, Inversor OLT e Macro
    double[] taxas_de_falha={2000,256,191,345,2381,120,50,500,2381,1075,667,191,381,2000}; //(fit)
    double[] tempo_de_reparo={6,4,6,6,22,4.25,4.25,22,2,2,6,6,8}; //(h)
    double[] tempo_de_inst={1,1,2.45,1.3,0,0.17,0.17,0,0.17,0.5,2.45,1.3,24}; //(h)
    double[] preco={0,0,379,600,0,SR*3.75,400,0,0,2700,379,20246,0}; //($)
    double[] KWh={70,0,0,0,0,0,0,0,0,0,0,800}; //(W)
    double[] areacobertura={0.25,0.1};//{micro,pico}//até 0.25
    double PkmVF=700; //instalar
    double PkmTF=57000; //trenching
    
    
    public SimuladorTCO_fotoMod()// throws IOException
    { 
    
    switch (typePon){
        case 0:
            preco[1]=220;
            KWh[1]=5;      
            preco[8]=3250;
            KWh[8]=16;
            break;
        case 1:
            preco[1]=930;
            KWh[1]=6.5;      
            preco[8]=13760;
            KWh[8]=34.4;
            break;
        case 2:
            preco[1]=1860;
            KWh[1]=13;      
            preco[8]=27510;
            KWh[8]=137.6;
            break;
        default:
            System.out.println("SR invalido");
    }
     
    //System.out.println(Arrays.toString(taxas_de_falha));
    //Numero de Equipamentos
    //{onu,dstep,spt,rn,fstep,oltc,fswt,pfstep,dswt,pdstep,GES,Micro,Pico,Ant,oltp,macro, olt splitter, olt splitter switch}
    //PICO, ONU, Painel ONU, Inversor ONU, Distribution Fiber Step, SPL, RN Chassi, Feeder Fiber Step, OLT, Painel OLT, Inversor OLT e Macro
    
    int [] numEq=new int[18];
    int numspq=0;
    int [][] CenarioMWCli;
    int [][] CenarioMWPEB;
    int [][] Cenaricli=CenarioBase();
    int clientemat=0,parclientemw=0;
    for(int i=0;i<(n*N);i++){
        for(int j=0;j<(n*N);j++){
            clientemat=0;
        }
    }
    //Calculo de uso dos equipamentos
    int [][] Resposta = Hexagono (areacobertura[1]);
    numEq[0] = Resposta.length;
    numEq[1] = Resposta.length;
    numEq[2] = 0;
    if (typeSF ==1){
        switch(typePon){
            case 0:
                if(SR ==8){
                numEq[2] = 76;
                }else {
                numEq[2] = 44;
                }
                break;
            case 1:
                if(SR ==8){
                numEq[2] = 162;
                }else {
                numEq[2] = 89;
                }
                break;
            case 2:
                if(SR ==8){
                numEq[2] = 647;
                }else {
                numEq[2] = 360;
                }
                break;
            default:
                System.out.println("SR invalido");
        }
    } else {
    numEq[2] = 0;
    }
    if (typeSF==1){
    numEq[3] = Resposta.length;
    } else {
    numEq[3] = 0;
    }
    numEq[4] = 0;
    
    int auxAtual=0;
    for (int i=0;i<Resposta.length;i++){
        if(auxAtual == Resposta[1][2]){
            numEq[4]++;
        } else {
            auxAtual++;
        }
    }
    
    numEq[5] = Resposta[Resposta.length-1][2];
    numEq[6] = Resposta[Resposta.length-1][2];
    numEq[7] = (int)Math.ceil((Math.sqrt(Resposta[Resposta.length-1][2])+1)*Math.sqrt(Resposta[Resposta.length-1][2]));
    numEq[8] = Resposta[Resposta.length-1][2]+1;
    numEq[9] = (int)Math.ceil(numEq[6]/NPC);
    numEq[10] = 0;//num de inversores temporarias
    if(typeSF==1){
        switch(typePon){
            case 0:
                if(SR ==8){
                numEq[10] = 76;
                }else {
                numEq[10] = 44;
                }
                break;
            case 1:
                if(SR ==8){
                numEq[10] = 162;
                }else {
                numEq[10] = 89;
                }
                break;
            case 2:
                if(SR ==8){
                numEq[10] = 647;
                }else {
                numEq[10] = 360;
                }
                break;
            default:
                System.out.println("SR invalido");
        }
    } else {
    numEq[10]=0;
    }
    
    numEq[11] = 1;//num de placas temporarias
            if (typeSF ==1){
                if (SR ==8 && typePon == 2){
                    numEq[11] = 2;
                }
            } else {
                numEq[11]=0;
            }
    numEq[12] = 1; 
    
    
           
    //Calculo do num total de equipamentos
    int numtoteq=0;
    for(int i=0;i<18;i++){
        numtoteq=numtoteq+numEq[i];
    }
    Equipamento [] Equi = new Equipamento[numtoteq];
    for(int i=0;i<numtoteq;i++){
        Equi[i] = new Equipamento();
    }
    int parpas=0;
    //Sem clientes e proteção, cenario 10x10 e distancia 1/24
    //Atualizar valores de numEq para cada equipamento
    //det distancias
    double [][] respostaSPL = new double [Resposta.length][2];
    int pc = (N*n-1)/2;
    auxAtual = 1;
    double [] auxResposta = new double [Resposta[Resposta.length-1][2]];
    int [][] posSPL = new int [Resposta[Resposta.length-1][2]][2];
    
    for (int i=0;i<Resposta.length;i++){
        respostaSPL[i][1] = Resposta[i][2];
        respostaSPL[i][2] = Math.sqrt((Resposta[i][0]-pc)*(Resposta[i][0]-pc)
                +(Resposta[i][1]-pc)*(Resposta[i][1]-pc));
    }
    for (int i=0;i<auxResposta.length;i++){
        auxResposta[i] = Integer.MAX_VALUE;
    }
    for(int i=0;i<respostaSPL.length;i++){
        if(auxAtual != respostaSPL[i][1]){
            auxAtual++;
        }
        if(respostaSPL[i][0]<auxResposta[auxAtual-1]){
            posSPL[auxAtual-1][0] = Resposta[i][0];
            posSPL[auxAtual-1][1] = Resposta[i][1];
            auxResposta[auxAtual-1]=respostaSPL[i][0];
        }
    }
    auxAtual=1;
    double [] DisestB= new double[Resposta.length];
    double [] DisestBali= new double[Resposta.length];
    double [] DisestBdist= new double[Resposta.length];
    double [] DisSPL= new double[Resposta[Resposta.length-1][2]];
    for (int i=0;i<Resposta.length;i++){
        if(auxAtual != respostaSPL[i][1]){
            auxAtual++;
        }
        DisestBdist[i]=Math.sqrt(Math.pow((Resposta[i][0]-posSPL[auxAtual-1][0]),2)+Math.pow((Resposta[i][1]-posSPL[auxAtual-1][1]),2));
        DisestBali[i]=Math.sqrt(Math.pow(posSPL[auxAtual-1][0]-pc,2)+Math.pow(posSPL[auxAtual-1][1]-pc,2));
        DisestB[i]=DisestBdist[i]+DisestBali[i];
    }
    //Pico
    if(numEq[0]>0){
        for(int i=0;i<numEq[0];i++){
            Equi[parpas].distancia= DisestB[i];
            System.out.println(Equi[parpas].distancia);
            Equi[parpas].taxadefalha=taxas_de_falha[0]/1000000000.0;
            Equi[parpas].tempodereparo=tempo_de_reparo[0];
            Equi[parpas].instalacao=tempo_de_inst[0];
            Equi[parpas].custoeq=preco[0];
            Equi[parpas].ConsumoEnergia=KWh[0];
            Equi[parpas].tipo=1;
            parpas++;
        }
    }
    //ONU
    if(numEq[1]>0){
        for(int i=0;i<numEq[1];i++){
            Equi[parpas].distancia= DisestB[i];
            Equi[parpas].taxadefalha=taxas_de_falha[1]/1000000000.0;
            Equi[parpas].tempodereparo=tempo_de_reparo[1];
            Equi[parpas].instalacao=tempo_de_inst[1]+ Equi[parpas].distancia/vel;
            Equi[parpas].custoeq=preco[1];
            Equi[parpas].ConsumoEnergia=KWh[1];
            Equi[parpas].tipo=2;
            parpas++;
        }
    }
    //PN ONU
    if(numEq[2]>0){
        for(int i=0;i<numEq[2]/2;i++){
            Equi[parpas].distancia= DisestB[i];
            Equi[parpas].taxadefalha=taxas_de_falha[2]/1000000000.0;
            Equi[parpas].tempodereparo=tempo_de_reparo[2];
            Equi[parpas].instalacao=tempo_de_inst[2] + Equi[parpas].distancia/vel;
            Equi[parpas].custoeq=preco[2];
            Equi[parpas].ConsumoEnergia=KWh[2];
            Equi[parpas].tipo=3;
            parpas++;
        }
        for(int i=0;i<numEq[2]/2;i++){
            Equi[parpas].distancia= DisestB[i];
            Equi[parpas].taxadefalha=taxas_de_falha[2]/1000000000.0;
            Equi[parpas].tempodereparo=tempo_de_reparo[2];
            Equi[parpas].instalacao=tempo_de_inst[2] + Equi[parpas].distancia/vel;
            Equi[parpas].custoeq=preco[2];
            Equi[parpas].ConsumoEnergia=KWh[2];
            Equi[parpas].tipo=3;
            parpas++;
        }
    }
    //InV ONU
    if(numEq[3]>0){
        for(int i=0;i<numEq[3];i++){
            Equi[parpas].distancia=DisestB[i];
            Equi[parpas].taxadefalha=taxas_de_falha[3]/1000000000.0;
            Equi[parpas].tempodereparo=tempo_de_reparo[3];
            Equi[parpas].instalacao=tempo_de_inst[3] + Equi[parpas].distancia/vel;
            Equi[parpas].custoeq=preco[3];
            Equi[parpas].ConsumoEnergia=KWh[3];
            Equi[parpas].tipo=4;
            parpas++;
        }
    }
    //Distribuition Fiber Step
    double medclusters= 0;
    for (int i=0;i<Resposta.length;i++){
        medclusters=medclusters+DisestB[i];
    }
    int medEstB = Resposta.length/Resposta[Resposta.length-1][2];
    double dH = 3*areacobertura[1];
    double dV = Math.sqrt(3)*areacobertura[1];
    double sqrtMedEstB = Math.sqrt(medEstB);
    double auxDiv = (medEstB - sqrtMedEstB) + sqrtMedEstB -1;
    double raioM = ((medEstB - sqrtMedEstB)*dH + (sqrtMedEstB-1)*dV)/auxDiv;
    
    medclusters=medclusters/numEq[0];
    
    double sumDisEst= 0;
        
    if(numEq[4]>0){
        for(int i=0;i<numEq[4];i++){
            Equi[parpas].distancia=medclusters;
            Equi[parpas].taxadefalha=raioM*areacobertura[1]*taxas_de_falha[4]/1000000000.0;
            Equi[parpas].tempodereparo=tempo_de_reparo[4];
            Equi[parpas].instalacao=tempo_de_inst[4];
            Equi[parpas].custoeq=preco[4];
            Equi[parpas].ConsumoEnergia=KWh[4];
            Equi[parpas].tipo=5;
            parpas++;
        }
    }
    
    //SPL
    if(numEq[5]>0){
        for(int i=0;i<numEq[5];i++){
            Equi[parpas].distancia=Math.sqrt((respostaSPL[i][0]-pc)*(respostaSPL[i][0]-pc)
                +(respostaSPL[i][1]-pc)*(respostaSPL[i][1]-pc));
            sumDisEst += Equi[parpas].distancia; 
            Equi[parpas].taxadefalha=taxas_de_falha[5]/1000000000.0;
            Equi[parpas].tempodereparo=tempo_de_reparo[5];
            Equi[parpas].instalacao=tempo_de_inst[5] + Equi[parpas].distancia/vel;
            Equi[parpas].custoeq=preco[5];
            Equi[parpas].ConsumoEnergia=KWh[5];
            Equi[parpas].tipo=6;
            parpas++;
        }
    }
    
    //RN chassi
    if(numEq[6]>0){
        for(int i=0;i<numEq[6];i++){
            Equi[parpas].distancia=Math.sqrt((respostaSPL[i][0]-pc)*(respostaSPL[i][0]-pc)
                +(respostaSPL[i][1]-pc)*(respostaSPL[i][1]-pc));
            Equi[parpas].taxadefalha=taxas_de_falha[6]/1000000000.0;
            Equi[parpas].tempodereparo=tempo_de_reparo[6];
            Equi[parpas].instalacao=tempo_de_inst[6] + Equi[parpas].distancia/vel;
            Equi[parpas].custoeq=preco[6];
            Equi[parpas].ConsumoEnergia=KWh[6];
            Equi[parpas].tipo=7;
            parpas++;
        }
    }
    
    //Feeder Fiber Step
    medEstB = Resposta[Resposta.length-1][2];
    
    switch (SR){
            case 16:
                dH = 10.5/2;
                dV = 7 * Math.sqrt(3)/4;
                break;
            case 32:
                dH = 10.5/2;
                dV = 7 * Math.sqrt(3)/2;
                break;
            case 64:
                dH = 10.5;
                dV = 7 * Math.sqrt(3)/2;
                break;
            default:
                break;
    }
    
    sqrtMedEstB = Math.sqrt(medEstB);
    auxDiv = (medEstB - sqrtMedEstB) + sqrtMedEstB -1;
    raioM = ((medEstB - sqrtMedEstB)*dH + (sqrtMedEstB-1)*dV)/auxDiv;
    medclusters = 0;
    
    for (int i=0; i<numEq[5]; i++){
        medclusters +=Math.sqrt((respostaSPL[i][0]-pc)*(respostaSPL[i][0]-pc)
                +(respostaSPL[i][1]-pc)*(respostaSPL[i][1]-pc)); 
    }
    medclusters = medclusters/numEq[5];
    
    if(numEq[7]>0){
        for(int i=0;i<numEq[7];i++){
            Equi[parpas].distancia=medclusters;
            Equi[parpas].taxadefalha=raioM*areacobertura[1]*taxas_de_falha[7]/1000000000.0;
            Equi[parpas].tempodereparo=tempo_de_reparo[7];
            Equi[parpas].instalacao=tempo_de_inst[7];
            Equi[parpas].custoeq=preco[7];
            Equi[parpas].ConsumoEnergia=KWh[7];
            Equi[parpas].tipo=8;
            parpas++;
        }
    }
    
    //OLT atualizar os valores
    if(numEq[8]>0){
        for(int i=0;i<numEq[8];i++){
            Equi[parpas].distancia=0;
            Equi[parpas].taxadefalha=taxas_de_falha[8]/1000000000.0;
            Equi[parpas].tempodereparo=tempo_de_reparo[8];
            Equi[parpas].instalacao=tempo_de_inst[8] + Equi[parpas].distancia/vel;
            Equi[parpas].custoeq=preco[8];
            Equi[parpas].ConsumoEnergia=KWh[8];
            Equi[parpas].tipo=9;
            parpas++;
        }
    }
    //OLT chassi
    if(numEq[9]>0){
        for(int i=0;i<numEq[9];i++){
            Equi[parpas].distancia=0;
            Equi[parpas].taxadefalha=taxas_de_falha[9]/1000000000.0;
            Equi[parpas].tempodereparo=tempo_de_reparo[9];
            Equi[parpas].instalacao=tempo_de_inst[9] + Equi[parpas].distancia/vel;
            Equi[parpas].custoeq=preco[9];
            Equi[parpas].ConsumoEnergia=KWh[9];
            Equi[parpas].tipo=10;
            parpas++;
        }
    }
    //PN OLT
    if(numEq[10]>0){
        for(int i=0;i<numEq[10];i++){
            Equi[parpas].distancia=0;
            Equi[parpas].taxadefalha=taxas_de_falha[10]/1000000000.0;
            Equi[parpas].tempodereparo=tempo_de_reparo[10];
            Equi[parpas].instalacao=tempo_de_inst[10] + Equi[parpas].distancia/vel;
            Equi[parpas].custoeq=preco[10];
            Equi[parpas].ConsumoEnergia=KWh[10];
            Equi[parpas].tipo=11;
            parpas++;
        }
    }
    //InV OLT
    if(numEq[11]>0){
        for(int i=0;i<numEq[11];i++){
            Equi[parpas].distancia=0;
            Equi[parpas].taxadefalha=taxas_de_falha[11]/1000000000.0;
            Equi[parpas].tempodereparo=tempo_de_reparo[11];
            Equi[parpas].instalacao=tempo_de_inst[11];
            Equi[parpas].custoeq=preco[11] + Equi[parpas].distancia/vel;
            Equi[parpas].ConsumoEnergia=KWh[11];
            Equi[parpas].tipo=12;
            parpas++;
        }
    }
    
    //Macro
    if(numEq[12]>0){
        Equi[parpas].taxadefalha=taxas_de_falha[12]/1000000000.0;
        Equi[parpas].tempodereparo=tempo_de_reparo[12];
        Equi[parpas].instalacao=tempo_de_inst[12];
        Equi[parpas].custoeq=preco[12];
        Equi[parpas].ConsumoEnergia=KWh[12];
        Equi[parpas].tipo=13;
        parpas++;
    }
    //Taxa de Falha Geral
        double TXFG=0;
        for(int i=0;i<numtoteq;i++){
            TXFG=TXFG+Equi[i].taxadefalha;
        }
        System.out.println("taxa de falha geral:  "+TXFG);
        //System.out.println(TXFG);
    //tempo e numero de iterações
        
    //numero de falhas de ONU , FDS e Olt Chassi em Ny
        double TEFOnu=numEq[1]*((365*24*Ny)/(1000000000.0/taxas_de_falha[0]));
        //System.out.println("Onu"+TEFOnu);
        // Montecarlo e Resolução do Modelo de cutos de markov
        
        Random gerador = new Random();
        double SCustoSLA=0;
        double SCustoEQT=0;
        double SqCustoSLA=0;
        double SqCustoEQT=0;
        double DPSomProdS=0;
        double DPSomProdE=0;
        double somanumit = 0;
        double[] custoits= new double[Num_Tentativas];
        double[] custoite= new double[Num_Tentativas];
        double custoaluguelexe=0;
        double[] Eq1=new double[18];
        double CustoEnergiadesp = 0;
        double CustoEnergiadespBack=0;
        double Custo_Energia=0;
        for(int i=0;i<Num_Tentativas;i++){
            int numit=0;
            double SomProdS=0;
            double SomProdE=0;
            double SomaAluguel=0;
            double u = gerador.nextDouble();
            int[] Vetor_Equi_fal=new int[numtoteq];
            int FimWhile=0;
            int[] Eq2=new int[18];
            while(FimWhile==0){
                //System.out.println(Arrays.toString(Vetor_Equi_fal));
                double Lambda;
                int Num_i=0;
                int Contfalha=0;
                int[] posfalha=new int[4];
                int possel = 0;
                int posus1 = 0;
                double ProdutoCEXT;
                double ProdutoCSXT;
                double Produtoaluguel = 0;
                double CustoAlugueletapa=0;
                for(int j=0;j<numtoteq;j++){
                    if(Vetor_Equi_fal[j]==1){
                        if(Num_i<4){
                        Num_i++;
                        }
                        posfalha[Contfalha]=j;
                        if(Contfalha<3){
                        Contfalha++;
                        }
                    }
                }
                if(Num_i==0){
                    ProdutoCEXT=0;
                    ProdutoCSXT=0;
                    Lambda=TXFG;
                }else{
                    Equipamento [] Equifal = new Equipamento[Num_i];
                     for(int j=0;j<Num_i;j++){
                     Equifal[j] = new Equipamento(); 
                    }
                    for(int j=0;j<Num_i;j++){
                    Equifal[j]=Equi[posfalha[j]];
                    }
                    for(int j=0;j<Num_i;j++){
                    if(Equi[posfalha[j]].tipo==0){
                        Eq2[0]++;
                    }else if(Equi[posfalha[j]].tipo==1){
                        Eq2[1]++;
                    }else if(Equi[posfalha[j]].tipo==2){
                        Eq2[2]++;
                    }else if(Equi[posfalha[j]].tipo==3){
                        Eq2[3]++;
                    }else if(Equi[posfalha[j]].tipo==4){
                        Eq2[4]++;
                    }else if(Equi[posfalha[j]].tipo==5){
                        Eq2[5]++;
                    }else if(Equi[posfalha[j]].tipo==6){
                        Eq2[6]++;
                    }else if(Equi[posfalha[j]].tipo==7){
                        Eq2[7]++;
                    }else if(Equi[posfalha[j]].tipo==8){
                        Eq2[8]++;
                    }else if(Equi[posfalha[j]].tipo==9){
                        Eq2[9]++;
                    }else if(Equi[posfalha[j]].tipo==10){
                        Eq2[10]++;
                    }else if(Equi[posfalha[j]].tipo==11){
                        Eq2[11]++;
                    }else if(Equi[posfalha[j]].tipo==12){
                        Eq2[12]++;
                    }else if(Equi[posfalha[j]].tipo==13){
                        Eq2[13]++;
                    }else if(Equi[posfalha[j]].tipo==14){
                        Eq2[14]++;
                    }else if(Equi[posfalha[j]].tipo==15){
                        Eq2[15]++;
                    }else if(Equi[posfalha[j]].tipo==16){
                        Eq2[16]++;
                    }else if(Equi[posfalha[j]].tipo==17){
                        Eq2[17]++;
                    }
                    
                    }
                    double[][] Matriz_Parametro= new double[9][4];
                    //System.out.println(Num_i);
                    //System.out.println(Arrays.toString(VetorcUSTO));
                    double[] VetorRep=new double[Num_i];
                    for(int j=0;j<Num_i;j++){
                        VetorRep[j]=Equifal[j].tempodereparo; 
                    }
                    possel=0;
                    double ini=0;
                    double sumtxfds=0;
                    for(int j=0;j<Num_i;j++){
                        sumtxfds=sumtxfds+Equi[posfalha[j]].taxadefalha;
                        if(ini<VetorRep[j]){
                            ini=VetorRep[j];
                            possel=j;
                        }
                        }
                    Lambda=TXFG+(1/(Equi[posfalha[possel]].tempodereparo+(Equi[posfalha[possel]].distancia/vel)))-sumtxfds;
                    double Tempo=(1/Lambda);//*(Math.log(u));
                    //não generico
                    double CustoET;
                    
                    int numfalha;
                    double Custoret=0;
                    Custoret=(Equi[posfalha[possel]].tempodereparo)*PorcRep;
                    CustoET=sal+Custoret;
                    ProdutoCEXT=CustoET*Tempo;
                    for(int j=0;j<Num_i;j++){
                        if(Equi[j].ConsumoEnergia>0){
                            CustoEnergiadesp=CustoEnergiadesp+((Equi[j].ConsumoEnergia/1000.0)*PrKw*Tempo);
                        }
                    }
                }
                
                //SomProd
                SomProdE=SomProdE+ProdutoCEXT;
                //Seleciona novo estado
                double seletor = gerador.nextDouble();
                double seletorsum = seletor*Lambda;
                double SSF=0;
                if(Num_i==0){
                    double SSFA=0;
                    for(int j=0;j<numtoteq;j++){
                        SSF=SSF+Equi[j].taxadefalha;
                        
                        if(seletorsum>SSFA){
                            if(seletorsum<=SSF){
                                Vetor_Equi_fal[j]=1;
                            }
                        }
                        SSFA=SSF;
                    }
                }else{
                    //System.out.println("seletor" + seletorsum);
                    double SSFA=0;
                    for(int j=0;j<numtoteq;j++){
                        if(j==posfalha[possel]){
                            SSF=SSF+(1/(Equi[j].tempodereparo+(Equi[j].distancia/vel)));
                            if(seletorsum>SSFA){
                                if(seletorsum<=SSF){
                                    Vetor_Equi_fal[j]=0;
                                }
                            }
                            SSFA=SSF;
                        }else{
                        int pfi=0;    
                        for(int k=0;k<Num_i;k++){
                            if(j==posfalha[k]&&Num_i<=4){
                                pfi=1;
                            }
                        }
                        if(pfi==0){
                            SSF=SSF+Equi[j].taxadefalha;
                             if(seletorsum>SSFA){
                                if(seletorsum<=SSF){
                                    Vetor_Equi_fal[j]=1;
                                }
                            }
                            SSFA=SSF; 
                        }else{
                            SSFA=SSF;
                        }
                        }
                    }
                }
                int Num_j=0;
                numit++;
                //System.out.println(Arrays.toString(Eq2));
                if(Eq2[2]>=TEFOnu){
                    FimWhile=1;
                }
                }
                for(int j=0;j<18;j++){
                    Eq1[j]=Eq2[j]+Eq1[j];
                }
                //System.out.println(i);
                somanumit=somanumit+numit;
                //SomProd
                custoaluguelexe=custoaluguelexe+SomaAluguel;
                custoits[i]=SomProdS;
                custoite[i]=SomProdE;
                System.out.println("iteração: "+i);
            }
        for(int j=0;j<16;j++){
                    Eq1[j]=Eq1[j]/Num_Tentativas;
                }
        double medianumit=somanumit/Num_Tentativas;
        double CReparo;
        double SomaEQT=0;
        for(int i=0;i<Num_Tentativas;i++){
            SomaEQT=SomaEQT+custoite[i];
        }
        //OPEX
        
        //Reparo
        CReparo=SomaEQT/Num_Tentativas;
        
        //Energia
        double CEnergia=0;
        for(int i=0;i<13;i++){
                CEnergia=numEq[i]*(KWh[i]/1000.0)*(24-(5.2*typeSF)*365)*Ny * PrKw;
        }
        CEnergia=CEnergia-(CustoEnergiadesp/Num_Tentativas);
        
        //System.out.println(Arrays.toString(custoits));
        //.println(Arrays.toString(Eq1));
        
        //Aluguel de espaco para as celulas
        
        
        cAluguel=(numEq[0]+numEq[2]+numEq[10]) * prAluguel;
        
        //Custo com gerenciamento com 4 equipes compostas por 1 eng e 2 tec
        cGerenc = 12.300*4*12 * Ny;
        
        
        
        
        // CAPEX
        
        //Custo de Equipamentos
        double[] VetorCustoCapexEqui=new double[12];
        
        for (int i=0;i<13;i++){
        VetorCustoCapexEqui[i]=numEq[i]*preco[i];
        }
        double CustoCompraEqui=0.0;
        for(int i=0;i<13;i++){
            CustoCompraEqui=CustoCompraEqui+VetorCustoCapexEqui[i];
        }
        // Custo de instalacao 
        double custoTotInst=0;
        
        for (int i=0;i<13;i++){
            custoTotInst = Equi[parpas].instalacao;
        }
        
        //Custos de Fibra    
        double cFibra = 0;
        double distFiberSR8 = 268.71;
        double distFiberSR16 = 280.74;
        for (int i =0; i < DisestBdist.length; i++){
        sumDisEst += DisestBdist[i];
        }
        
        if (SR == 8){
            cFibra = (sumDisEst * PkmVF) + (distFiberSR8 * PkmTF);
        } else if (SR == 16){
            cFibra = (sumDisEst * PkmVF) + (distFiberSR16 * PkmTF);
        } else {
            System.out.println("SR nao ajustado");
        }
        
        ///Resultado
        String caract1="";
        String caract2="";
        String tiprot3="";
        String tiprot4="";
        String tiprot5="";
        String timw="";
        switch (typePon){
            case 0:
                caract1 = "GPON";
                break;
            case 1:
                caract1 = "10GPON";
                break;
            case 2:
                caract1 = "40GPON";
                break;
            default:
                
        }
        switch (SR){
            case 8:
                caract2 = "SR = 8";
                break;
            case 16:
                caract2 = "SR = 16";
                break;
            case 32:
                caract2 = "SR = 32";
                break;
            case 64:
                caract2 = "SR = 64";
                break;
            default:
        }
        ///Resultado
        //Ordem dos Resultados: CAPEX |Instalacao - Equipamento - Fibra (CAPEX)|  OPEX  | Reparo - Energia - Aluguel por m2 - Aluguel de Gerenciamento |
        System.out.println( caract1 +" , "+ caract2 + " : \n");
        System.out.println(CustoCompraEqui+custoTotInst+cFibra);                
        System.out.println(CustoCompraEqui);
        System.out.println(custoTotInst);
        System.out.println(cFibra);
        System.out.println(CEnergia+CReparo+cAluguel+cGerenc);
        System.out.println(CReparo);
        System.out.println(CEnergia);
        System.out.println(cAluguel);
        System.out.println(cGerenc);
        System.out.println( "Média de falhas de equipamentos:  "+"{PICO, ONU, Painel ONU, Inversor ONU, Distribution Fiber Step, SPL, RN Chassi, Feeder Fiber Step, OLT, OltChassi, Painel OLT, Inversor OLT e Macro}");
        System.out.println( "Média de falhas de equipamentos:  " + Arrays.toString(Eq1));
        
//escrever novo metodo de salvar
        
        
  

        
    }
    private class Equipamento{
    int tipo;
    double distancia;
    int[] clientes=new int[2];
    int[] protecao=new int[3];
    double taxadefalha;
    double tempodereparo;
    double custoeq;
    double instalacao;
    double ConsumoEnergia;
    }
    private double[] CriarONU(int i, int paronu){
        double[] Retorno={0,0,0,0,0,0,0,0,0,0};
        //preco,Qaux,PVaux,PHaux,Cli,Cli,Dis,TXf,TemR,TemI.
        int Quadraaux=(int) Math.floor((double) i/Math.pow(n,2));
        int PosAuxFFV=(int) Math.floor((double) Quadraaux/N);
        int PosAuxFFH= Quadraaux-N*PosAuxFFV;
        int Posvetaux=(int) Math.floor((double) paronu/n);
        int Poshoraux= paronu-n*Posvetaux;
        double distanciadaux=0.0; //distancia DL
        double distanciafaux=0.0; //distancia FL
        Retorno[0]=preco[0];
        Retorno[1]=Quadraaux;
        Retorno[2]=Posvetaux;
        Retorno[3]=Poshoraux;  
        Retorno[4]=i; 
        Retorno[5]=i;
        if(PosAuxFFV<N/2.0){
            if(PosAuxFFH<N/2.0){
                distanciafaux=(((N-(2*PosAuxFFH+1))+(N-(2*PosAuxFFV+1)))*(L/2.0));
            }else{
                distanciafaux=(((N-(2*(N-PosAuxFFH-1)+1))+(N-(2*PosAuxFFV+1)))*(L/2.0));
            }
        }else{
            PosAuxFFV=N-PosAuxFFV-1;
            if(PosAuxFFH<N/2.0){
                distanciafaux=(((N-(2*PosAuxFFH+1))+(N-(2*PosAuxFFV+1)))*(L/2.0));
            }else{
                distanciafaux=(((N-(2*(N-PosAuxFFH-1)+1))+(N-(2*PosAuxFFV+1)))*(L/2.0));
            }
        }
        if(Posvetaux<n/2){
            if(Poshoraux<n/2){
                distanciadaux=(((n-(2*Poshoraux+1))+(n-(2*Posvetaux+1)))*(l/2.0));
            }else{
                distanciadaux=(((n-(2*(n-Poshoraux-1)+1))+(n-(2*Posvetaux+1)))*(l/2.0));
            }
        }else{
            Posvetaux=n-Posvetaux-1;
            if(Poshoraux<n/2){
                distanciadaux=(((n-(2*Poshoraux+1))+(n-(2*Posvetaux+1)))*(l/2.0));
            }else{
                distanciadaux=(((n-(2*(n-Poshoraux-1)+1))+(n-(2*Posvetaux+1)))*(l/2.0));
            }
        }
        Retorno[6]=(distanciadaux+distanciafaux);
        Retorno[7]=taxas_de_falha[0]/1000000000.0;
        Retorno[8]=tempo_de_reparo[0];  
        Retorno[9]=tempo_de_inst[0];
        
        return Retorno;
    }
    private double[] CriarDisStep(int i,int parpas, int pards, int numcli,int[][]Onucli){
        double[] Retorno={0,0,0,0,0,0,0,0};
        //Cli,Cli,Dis,Prot,Prot,Prot,TXf,TemR.
        int Quadraaux=(int) Math.floor((double) i/((n+1)*n));
        int PosAuxFFV=(int) Math.floor((double) Quadraaux/N);
        int PosAuxFFH= Quadraaux-N*PosAuxFFV;
        int Posvetaux=(int) Math.floor((double) (pards/(n+1)));
        int Posvetaux2=Posvetaux;
        int Poshoraux= pards-(n+1)*Posvetaux;
        double distanciadaux=0.0;
        double distanciafaux=0.0;
        Retorno[3]=4;
        if(PosAuxFFV<N/2.0){
            if(PosAuxFFH<N/2.0){
                distanciafaux=(((N-(2*PosAuxFFH+1))+(N-(2*PosAuxFFV+1)))*(L/2.0));
            }else{
                distanciafaux=(((N-(2*(N-PosAuxFFH-1)+1))+(N-(2*PosAuxFFV+1)))*(L/2.0));
            }
        }else{
            PosAuxFFV=N-PosAuxFFV-1;
            if(PosAuxFFH<N/2.0){
                distanciafaux=(((N-(2*PosAuxFFH+1))+(N-(2*PosAuxFFV+1)))*(L/2.0));
            }else{
                distanciafaux=(((N-(2*(N-PosAuxFFH-1)+1))+(N-(2*PosAuxFFV+1)))*(L/2.0));
            }
        }
        if(Posvetaux<(n/2)){
            if(Poshoraux<((n-2)/2)){
                distanciadaux=((n-(2*Posvetaux+1))*(l/2.0))+((n/2)-1-Poshoraux)*l;
            }else if(Poshoraux>((n+2)/2)){
                distanciadaux=((n-(2*Posvetaux+1))*(l/2.0))+(Poshoraux-(n/2)-1)*l;
            }else{
                distanciadaux=(n-(2*Posvetaux+1))*(l/2.0);
            }
        }else{
            if(Poshoraux<((n-2)/2)){
                distanciadaux=((n-(2*(n-Posvetaux-1)+1))*(l/2.0))+((n/2)-1-Poshoraux)*l;
            }else if(Poshoraux>((n+2)/2)){
                distanciadaux=((n-(2*(n-Posvetaux-1)+1))*(l/2.0))+(Poshoraux-(n/2)-1)*l;
            }else{
                distanciadaux=(n-(2*(n-Posvetaux-1)+1))*(l/2.0);
            }
        }
        Retorno[2]=(distanciadaux+distanciafaux);
        int[] VetorClDF = new int[numcli];
        for(int j=0;j<numcli;j++){
            if(Quadraaux==Onucli[0][j]){
            if(Poshoraux==n/2){
                if(Posvetaux<n/2){
                    if(Onucli[1][j]<=Posvetaux){
                        VetorClDF[j]=1;
                    }
                }else{
                    if(Onucli[1][j]>=Posvetaux){
                        VetorClDF[j]=1;
                    }
                }
            }else if(Poshoraux<n/2){
                if(Onucli[1][j]==Posvetaux){
                    if(Onucli[2][j]<=Poshoraux){
                        VetorClDF[j]=1;
                    }
                }   
            }else if(Poshoraux>n/2){
                if(Onucli[1][j]==Posvetaux){
                    if(Onucli[2][j]>=((Poshoraux)-1)){
                        VetorClDF[j]=1;
                    }
                }
            }
            }
        }
        int[] vetorprot=new int[2];
        int Q=0;
        if(Posvetaux2==((n/2.0)-1)){
            Retorno[3]=5;
            Q=(int) (Quadraaux*Math.pow(n,2));
            if(Poshoraux==(n/2.0)){
                vetorprot[0]=((n/2)*n);
                vetorprot[1]=((n*n)-1);
            }else if(Poshoraux<(n/2.0)){
                vetorprot[0]=((n/2)*n);
                vetorprot[1]=(Poshoraux)*(n/2)+((n/2)*n)+((n/2)-1);
            }else{
                vetorprot[0]=(Poshoraux-1)*(n/2)+((n/2)*n);
                vetorprot[1]=((n*n)-1);
            }
        }else if(Posvetaux2==(n/2.0)){
            Retorno[3]=5;
            Q=(int) (Quadraaux*Math.pow(n,2));
            if(Poshoraux==(n/2.0)){
                vetorprot[0]=0;
                vetorprot[1]=(((n/2)*n)-1);
            }else if(Poshoraux<(n/2.0)){
                vetorprot[0]=0;
                vetorprot[1]=Poshoraux*(n/2)+((n/2)-1);
            }else{
                vetorprot[0]=(Poshoraux-1)*(n/2);
                vetorprot[1]=(((n/2)*n)-1);
            }
        }
        Retorno[4]=vetorprot[0]+Q; 
        Retorno[5]=vetorprot[1]+Q;
        int[] clientes=CriarVCli(VetorClDF);
        Retorno[0]=clientes[0];
        Retorno[1]=clientes[1];
        Retorno[6]=(taxas_de_falha[1]*l)/1000000000.0;
        Retorno[7]=tempo_de_reparo[1];
        parpas++;
        return Retorno;
    }
    private double[] Criar_Splitter(int i, int parpas, int numspq, int numcli, int Cps, int CliQ){
    double[] Retorno={0,0,0,0,0,0,0,0,0};
    //Cli,Cli,Dis,preco,TXf,TemR,TemI,Cps.
    int Quadraaux=(int) Math.floor((double) i/numspq);
    int PosAuxFFV=(int) Math.floor((double) Quadraaux/N);
    int PosAuxFFH= Quadraaux-N*PosAuxFFV;
    double distanciafaux=0.0;
        if(PosAuxFFV<N/2.0){
            if(PosAuxFFH<N/2.0){
                distanciafaux=(((N-(2*PosAuxFFH+1))+(N-(2*PosAuxFFV+1)))*(L/2.0));
            }else{
                PosAuxFFH=N-PosAuxFFH-1;
                distanciafaux=(((N-(2*PosAuxFFH+1))+(N-(2*PosAuxFFV+1)))*(L/2.0));
            }
        }else{
            PosAuxFFV=N-PosAuxFFV-1;
            if(PosAuxFFH<N/2.0){
                distanciafaux=(((N-(2*PosAuxFFH+1))+(N-(2*PosAuxFFV+1)))*(L/2.0));
            }else{
                PosAuxFFH=N-PosAuxFFH-1;
                distanciafaux=(((N-(2*PosAuxFFH+1))+(N-(2*PosAuxFFV+1)))*(L/2.0));
            }
        }
        Retorno[2]=distanciafaux;
        int[] VetorClSP = new int[numcli];
        if(CliQ>SR){
            for(int j=Cps;j<(Cps+SR);j++){
                VetorClSP[j]=1;
            }
            CliQ=CliQ-SR;
            Cps=Cps+SR;
        }else{
            for(int j=Cps;j<(Cps+CliQ);j++){
                VetorClSP[j]=1;
            }
            Cps=Cps+CliQ;
            CliQ=(int) Math.pow(n,2);
        }
        int[] clientes=CriarVCli(VetorClSP);
        Retorno[0]=clientes[0];
        Retorno[1]=clientes[1];
        Retorno[4]=taxas_de_falha[2]/1000000000.0;
        Retorno[5]=tempo_de_reparo[2];
        Retorno[6]=tempo_de_inst[2];
        Retorno[3]=preco[2];
        Retorno[7]=Cps;
        Retorno[8]=CliQ;
        parpas++;
        return Retorno;
    }
    private double[] CriarRn_Chassi(int i,int CpRn,int RCliQ,int numcli){
        double[] Retorno={0,0,0,0,0,0,0,0};
        //Cli,Cli,Dis,preco,TXf,TemR,TemI,Cps.
        int Quadraaux=i;
        int PosAuxFFV=(int) Math.floor((double) Quadraaux/N);
        int PosAuxFFH= Quadraaux-N*PosAuxFFV;
        double distanciafaux=0.0;
        if(PosAuxFFV<N/2.0){
            if(PosAuxFFH<N/2.0){
                distanciafaux=(((N-(2*PosAuxFFH+1))+(N-(2*PosAuxFFV+1)))*(L/2.0));
            }else{
                PosAuxFFH=N-PosAuxFFH-1;
                distanciafaux=(((N-(2*PosAuxFFH+1))+(N-(2*PosAuxFFV+1)))*(L/2.0));
            }
        }else{
            PosAuxFFV=N-PosAuxFFV-1;
            if(PosAuxFFH<N/2.0){
                distanciafaux=(((N-(2*PosAuxFFH+1))+(N-(2*PosAuxFFV+1)))*(L/2.0));
            }else{
                PosAuxFFH=N-PosAuxFFH-1;
                distanciafaux=(((N-(2*PosAuxFFH+1))+(N-(2*PosAuxFFV+1)))*(L/2.0));
            }
        }
        Retorno[2]=distanciafaux;
        int[] VetorClRn = new int[numcli];
        for(int j=CpRn;j<(CpRn+RCliQ);j++){
                VetorClRn[j]=1;
        }
        CpRn=CpRn+RCliQ;
        int[] clientes=CriarVCli(VetorClRn);
        Retorno[0]=clientes[0];
        Retorno[1]=clientes[1];
        Retorno[4]=taxas_de_falha[4]/1000000000.0;
        Retorno[5]=tempo_de_reparo[4];
        Retorno[6]=tempo_de_inst[4];
        Retorno[3]=preco[4];
        Retorno[7]=CpRn;
        return Retorno;
    }
    private double[] CriarFedStep(int i,int parpas,int[][]Quadracli,int numcli,int NQ){
        double[] Retorno={0,0,0,0,0,0,0,0};
        //Cli,Cli,Dis,Prot,Prot,Prot,TXf,TemR.
        int Posvetaux=(int) Math.floor((double) (i/(N+1)));
        int Posvetaux2=Posvetaux;
        int Poshoraux= i-(N+1)*Posvetaux;
        double distanciafaux=0.0;
        Retorno[3]=1;
        if(Posvetaux<N/2){
            if(Poshoraux==N/2){
                distanciafaux=((N-2)/2-Posvetaux)*L;
            }else if(Poshoraux<N/2){
                distanciafaux=((N-(2*Posvetaux+1))*(L/2.0))+((N-2)/2-Poshoraux)*L;
            }else{
                distanciafaux=((N-(2*Posvetaux+1))*(L/2.0))+((N-2)/2-(N-Poshoraux))*L;
            }
        }else{
            Posvetaux=N-Posvetaux-1;
            if(Poshoraux==N/2){
                distanciafaux=((N-2)/2-Posvetaux)*L;
            }else if(Poshoraux<N/2){
                distanciafaux=((N-(2*Posvetaux+1))*(L/2.0))+((N-2)/2-Poshoraux)*L;
            }else{
                distanciafaux=((N-(2*Posvetaux+1))*(L/2.0))+((N-2)/2-(N-Poshoraux))*L;
            }
        }
        Retorno[2]=distanciafaux;
        int[] VetorClFF = new int[NQ];
        int[] VetorClFF2 = new int[numcli];
        int NCQ=(int) Math.pow(n,2);
        int contff=0;
        for (int j=0;j<NQ;j++){
            if(Poshoraux==N/2){
                if(Posvetaux<N/2){
                    if(Quadracli[0][j]<=Posvetaux){
                        VetorClFF[j]=1;
                    }
                }else{
                    if(Quadracli[0][j]>=Posvetaux){
                        VetorClFF[j]=1;
                    }
                }
            }else if(Poshoraux<N/2){
                if(Quadracli[0][j]==Posvetaux){
                    if(Quadracli[1][j]<=Poshoraux){
                        VetorClFF[j]=1;
                    }
                }   
            }else if(Poshoraux>N/2){
                if(Quadracli[0][j]==Posvetaux){
                if(Quadracli[1][j]>=((Poshoraux)-1)){
                        VetorClFF[j]=1;
                }
                }
            }
            for(int k=contff;k<(contff+NCQ);k++){
                if(VetorClFF[j]==1){
                    VetorClFF2[k]=1;
                }else{
                    VetorClFF2[k]=0;
                }
            }
            contff=contff+NCQ;
        }
        int[] vetorprot=new int[2];
        int Q=0;
        if(Posvetaux2==((N/2.0)-1)){
            Retorno[3]=2;
            if(Poshoraux==(n/2.0)){
                vetorprot[0]=((N/2)*N);
                vetorprot[1]=((N*N)-1);
            }else if(Poshoraux<(N/2.0)){
                vetorprot[0]=((N/2)*N);
                vetorprot[1]=(Poshoraux)*(N/2)+((N/2)*N)+((N/2)-1);
            }else{
                vetorprot[0]=(Poshoraux-1)*(N/2)+((N/2)*N);
                vetorprot[1]=((N*N)-1);
            }
        }else if(Posvetaux2==(N/2.0)){
            Retorno[3]=2;
            if(Poshoraux==(N/2.0)){
                vetorprot[0]=0;
                vetorprot[1]=(((N/2)*N)-1);
            }else if(Poshoraux<(N/2.0)){
                vetorprot[0]=0;
                vetorprot[1]=Poshoraux*(N/2)+((N/2)-1);
            }else{
                vetorprot[0]=(Poshoraux-1)*(N/2);
                vetorprot[1]=(((N/2)*N)-1);
            }
        }
        int CQ=n*n;
        Retorno[4]=vetorprot[0]*CQ;
        Retorno[5]=((vetorprot[1]+1)*CQ)-1;
        int[] clientes=CriarVCli(VetorClFF2);
        Retorno[0]=clientes[0];
        Retorno[1]=clientes[1];
        Retorno[6]=(taxas_de_falha[1]*L)/1000000000.0;
        Retorno[7]=tempo_de_reparo[1];
        return Retorno;
    }
    private double[] CriarOltC(int i,int alnumspt,int numcli,int numdissteps,int numspt,int contoc,Equipamento[] Equi){
        double[] Retorno={0,0,0,0,0,0,0};
        //Cli,Cli,Preco,TXf,TemR,TemI,contoc.
        int limoc;
        if(alnumspt<=NPC){
            limoc=alnumspt;
        }else{
            limoc=NPC;
            alnumspt=alnumspt-NPC;
        }  
        int[] VetorClOC = new int[numcli];
        for(int j=0;j<numcli;j++){
            VetorClOC[j]=0;
        }
        
        for(int j=contoc;j<contoc+limoc;j++){
            /*
            for(int k=0;k<numcli;k++){VetorClOCaux[k]=0;}
            for(int k=OLTPcli[j].clientes[0];k<=OLTPcli[j].clientes[1];k++){ 
                VetorClOCaux[k]=1;
            }*/
            VetorClOC=SomaVetores(ObterVetorCompleto(Equi[j+(numcli+numdissteps+numspt)].clientes), VetorClOC);
            
            
        }
        int[] clientes=CriarVCli(VetorClOC);
        contoc=contoc+limoc;
        Retorno[0]=clientes[0];
        Retorno[1]=clientes[1];
        Retorno[2]=preco[5];
        Retorno[3]=taxas_de_falha[5]/1000000000.0;
        Retorno[4]=tempo_de_reparo[5];
        Retorno[5]=tempo_de_inst[5];
        Retorno[6]=contoc;
        
        return Retorno;
    }
    private double[] CriarPFedStep(int i){
        double[] Retorno={0,0,0,0,0,0};
        //prot,prot,prot,dis,TXf,TemR.
        int Poshoraux=(int) Math.floor((double) (i/(N-1)));
        int Posvetaux= i-(N-1)*Poshoraux;
        double distanciafaux=0.0;
        if(Posvetaux==((N/2.0)-1)){
            if(Poshoraux<(N/2.0)){
                distanciafaux=(((N-1)/2.0)-Poshoraux)*L;
            }else{
                distanciafaux=(Poshoraux-((N-1)/2.0))*L;
            }
        }else if(Posvetaux>((N/2.0)-1)){
            if(Poshoraux<(N/2)){
                distanciafaux=(Posvetaux-((N-2)/2.0))*L+(((N-1)/2.0)-Poshoraux)*L;
            }else{
                distanciafaux=(Posvetaux-((N-2)/2.0))*L+(Poshoraux-((N-1)/2.0))*L;
            }
        }else{
            if(Poshoraux<(N/2.0)){
                distanciafaux=(((N-2)/2.0)-Posvetaux)*L+(((N-1)/2.0)-Poshoraux)*L;
            }else{
                distanciafaux=(((N-2)/2.0)-Posvetaux)*L+(Poshoraux-((N-1)/2.0))*L;
            }
        }
        Retorno[3]=distanciafaux;
        Retorno[0]=3;
        int Novo_N=n*(n);
        if(Posvetaux==((N/2.0)-1)){
            Retorno[1]=Poshoraux*(N*Novo_N);
            Retorno[2]=Poshoraux*(N*Novo_N)+(N*Novo_N-1);        
        }else if(Posvetaux<((N/2.0)-1)){
            Retorno[1]=Poshoraux*(N*Novo_N);
            Retorno[2]=Poshoraux*(N*Novo_N)+(Posvetaux*Novo_N+Novo_N-1);  
        }else{
            Retorno[1]=Poshoraux*(N*Novo_N)+(Posvetaux*Novo_N+Novo_N);
            Retorno[2]=Poshoraux*(N*Novo_N)+(N*Novo_N-1);    
        }
        Retorno[4]=(taxas_de_falha[1]*L)/1000000000.0;
        Retorno[5]=tempo_de_reparo[1];
        return Retorno;
    } 
    private double[] CriarPdisStep(int i,int pards){
        double[] Retorno={0,0,0,0,0,0};
        //prot,prot,prot,dis,TXf,TemR.
        int detquadra=(int) Math.floor((i/((n-1)*n)));
        int PosAuxFFH=(int) Math.floor((double) detquadra/N);
        int PosAuxFFV=detquadra-N*PosAuxFFH;
        int Poshoraux=(int) Math.floor((double) (pards/(n-1)));
        int Posvetaux= pards-(n-1)*Poshoraux;
        //determinar distancia
        double distanciaFL=0;
        double distanciaDL=0;
        double distancia=0;
        if(PosAuxFFV<N/2.0){
            if(PosAuxFFH<N/2.0){
                distanciaFL=(((N-(2*PosAuxFFH+1))+(N-(2*PosAuxFFV+1)))*(L/2.0));
            }else{
                PosAuxFFH=N-PosAuxFFH-1;
                distanciaFL=(((N-(2*PosAuxFFH+1))+(N-(2*PosAuxFFV+1)))*(L/2.0));
                }
            }else{
                PosAuxFFV=N-PosAuxFFV-1;
            if(PosAuxFFH<N/2.0){
                distanciaFL=(((N-(2*PosAuxFFH+1))+(N-(2*PosAuxFFV+1)))*(L/2.0));
            }else{
                PosAuxFFH=N-PosAuxFFH-1;
                distanciaFL=(((N-(2*PosAuxFFH+1))+(N-(2*PosAuxFFV+1)))*(L/2.0));
                }
            }
        if(Posvetaux==((n/2.0)-1)){
            if(Poshoraux<(n/2.0)){
                distanciaDL=(((n-1)/2.0)-Poshoraux)*l;
            }else{
                distanciaDL=(Poshoraux-((n-1)/2.0))*l;
            }
        }else if(Posvetaux>((n/2.0)-1)){
            if(Poshoraux<(n/2)){
                distanciaDL=(Posvetaux-((n-2)/2.0))*l+(((n-1)/2.0)-Poshoraux)*l;
            }else{
                distanciaDL=(Posvetaux-((n-2)/2.0))*l+(Poshoraux-((n-1)/2.0))*l;
            }
        }else{
            if(Poshoraux<(n/2.0)){
                distanciaDL=(((n-2)/2.0)-Posvetaux)*l+(((n-1)/2.0)-Poshoraux)*l;
            }else{
                distanciaDL=(((n-2)/2.0)-Posvetaux)*l+(Poshoraux-((n-1)/2.0))*l;
            }
        }
        distancia=distanciaFL+distanciaDL;
        Retorno[3]=distancia;
        //determinar clientes
        int[]ClientesP=new int[2];
        Retorno[0]=6;
        if(Posvetaux==((n/2.0)-1)){
            ClientesP[0]=Poshoraux*n;
            ClientesP[1]=Poshoraux*n+(n-1);        
        }else if(Posvetaux<((n/2.0)-1)){
            ClientesP[0]=Poshoraux*n;
            ClientesP[1]=Poshoraux*n+(Posvetaux);
        }else{
            ClientesP[0]=Poshoraux*n+(Posvetaux+1);
            ClientesP[1]=Poshoraux*n+(n-1);
        }
        int inicio=detquadra*(n*(n));
        Retorno[1]=inicio+ClientesP[0];
        Retorno[2]=inicio+ClientesP[1];
        Retorno[4]=(taxas_de_falha[1]*l)/1000000000.0;
        Retorno[5]=tempo_de_reparo[1];
        return Retorno;
    }
    private int[] SomaVetores(int []  va, int [] vb){
        int comp=(int) (n*n*N*N);
        int[] vc=new int[comp];
        for (int i=0; i< comp; i++ ){
            vc[i] = va[i] + vb[i];
        }
        return vc;
    }
    private int[] CriarVCli(int[] va){
        int[] vetau1 =new int[2];
        int numcli=(int) (Math.pow(n,2)*Math.pow(N,2));
        for(int j=1;j<numcli;j++){
            if(va[0]==1){
                if(va[(numcli-1)]==1){
                    vetau1[0]=0;
                    vetau1[1]=numcli-1;
                }else{
                    vetau1[0]=0;
                        if(va[j]==0){
                            if(va[j-1]==1){vetau1[1]=(j-1);}
                        }
            }
            }else{
                if(va[j]==1){
                    if(va[j-1]==0){vetau1[0]=j;}
                }
                if(va[j]==1){
                    if(va[(numcli-1)]==1){
                        vetau1[1]=(numcli-1);
                    }else{
                        if(va[j+1]==0){vetau1[1]=j;}
                    }
                }
            }
        }
        return vetau1;
    }
    private int[] ObterVetorCompleto(int[] va){
        int numcli=(int) (Math.pow(n,2)*Math.pow(N,2));
        int[] vb=new int[numcli];
        for(int k=0;k<numcli;k++){vb[k]=0;}
            
            for(int k=va[0];k<=va[1];k++){ 
                vb[k]=1;
            }
        return vb;
    }
    private int[] ObterVetorPrt2(int numcli,int[] vetorprot){
        int[] novo_vetor= new int[numcli];
        int nova_posicao;
        int Quadra;
        int Quadrante;
        int[] Quadranteaux=new int[2];
        int Pos1;
        int Pos2;
        //System.out.println("vetor prot" + Arrays.toString(vetorprot));
        int[] Pos2aux=new int[2];
        int[] Pos3aux=new int[2];
        int[] posfinal=new int[2];
        int n2=(int) Math.pow(n,2);
        for(int i=0;i<numcli;i++){
            Quadra=(int) Math.floor(i/n2);
            Pos1=i-(n2*Quadra);
            Quadrante=(int) Math.floor(Pos1/(n2/4.0));
            Quadranteaux[1]=(int) Math.floor(Quadrante/2.0);
            Quadranteaux[0]=Quadrante-Quadranteaux[1]*2;
            Pos2=(Pos1-(n2/4)*Quadrante);
            Pos2aux[0]=(int) Math.floor(Pos2/(n/2.0));;
            Pos2aux[1]=Pos2-Pos2aux[0]*(n/2);
            Pos3aux[0]=Pos2aux[0]+Quadranteaux[0]*(n/2);
            Pos3aux[1]=Pos2aux[1]+Quadranteaux[1]*(n/2);
            posfinal[0]=Pos3aux[1]*n+Pos3aux[0];
            posfinal[1]=posfinal[0]+(n2*Quadra);
            novo_vetor[posfinal[1]]=vetorprot[i];
        }
        //System.out.println(Arrays.toString(novo_vetor));
        return novo_vetor;
    }
    private int[] ObterVetorPrt(int numcli,int[] vetorprot){
        int[] novo_vetor= new int[numcli];
        int nova_posicao;
        int[] Quadra=new int[2];
        int[] Posaux=new int[2];
        int n2=(int) Math.pow(n,2);
        for(int i=0;i<numcli;i++){
            Quadra[0]=(int) Math.floor(i/n2);
            Quadra[1]=(int) ((Quadra[0]-(Math.floor(Quadra[0]/N)*N))*N+Math.floor(Quadra[0]/N));        
            Posaux[0]=i-(n2*Quadra[0]);
            Posaux[1]=(int) ((Posaux[0]-(Math.floor(Posaux[0]/n)*n))*n+Math.floor(Posaux[0]/n));
            nova_posicao=Posaux[1]+(n2*Quadra[1]);
            novo_vetor[nova_posicao]=vetorprot[i];
        }
        return novo_vetor;
    }
    private int[][] MatrizEstBase(double raio,int parmetodo){
        double EstBdis,EstBdismeioL,aux1,aux2=l+1;
        double NumeroEsp=0;
        EstBdis= 2*raio*Math.cos((Math.PI/6));
        EstBdismeioL=EstBdis*Math.cos((Math.PI/4));
        for(int i=0;i<N*n;i++){
            aux1=Math.abs((EstBdismeioL-(l*(i+1))));
            if(aux1==0){
                NumeroEsp=i+1;
            }else if(aux1<aux2){
                NumeroEsp=i+1;
            }
            aux2=aux1;
        }
        double DisUniEst= Math.ceil((n*N/NumeroEsp));
        if(DisUniEst % 2==0){
            DisUniEst=DisUniEst-1;
            //System.out.println("veio");
        }
        int NumEstB=(int) (Math.pow(Math.ceil(DisUniEst/2.0),2)+Math.pow(Math.floor(DisUniEst/2.0),2));
        //System.out.println("numero de estações base:"+NumEstB);
        int[][]cenario=new int[(n*N)][(n*N)];
        int NEBI=(int) Math.ceil(DisUniEst/2.0);
        int NEBP=(int) Math.floor(DisUniEst/2.0);
        int posini=1,contador1=1;
        
        //System.out.println(NumeroEsp);
        while(posini>=0){
            posini=(int) (((n*N)/2.0)-1-contador1*NumeroEsp);
            contador1++;
        }
        posini=(int) (posini+NumeroEsp);
        int poslin=posini;
        int poscol=posini;
        int parnp=0;
        for(int i=0;i<NumEstB;i++){
            cenario[poslin][poscol]=1;
            if((poscol+2*NumeroEsp<n*N)){
                poscol=(int) (poscol+2*NumeroEsp);
            }else{
                //i=i-1;
                if(parnp==0){
                    parnp=1;
                }else{
                    parnp=0;
                }
                poscol=(int) (posini+parnp*NumeroEsp);
                if((poslin+NumeroEsp<n*N)){
                    poslin=(int) (poslin+NumeroEsp);
                }
            }
        }
        int contador=0;
        int [][]PosEstB=new int[NumEstB][2];
        for(int i=0;i<n*N;i++){
            for(int j=0;j<n*N;j++){
                if(cenario[i][j]==1){
                    cenario[i][j]=contador;
                    PosEstB[contador][0]=i;
                    PosEstB[contador][1]=j;
                    contador++;
                }
            }
        }
        if(parmetodo==0){
        double[] DisPtEstB=new double[NumEstB];
        int[] SelEstB=new int[NumEstB];
        for(int i=0;i<n*N;i++){
            for(int j=0;j<n*N;j++){
                for(int k=0;k<NumEstB;k++){
                    DisPtEstB[k]=Math.sqrt(Math.pow((PosEstB[k][1]-j),2)+Math.pow((PosEstB[k][0]-i),2));
                }
                
                double Menordis=DisPtEstB[0];
                for(int k=0;k<NumEstB;k++){
                    if(DisPtEstB[k]<Menordis){
                        Menordis=DisPtEstB[k];
                    }
                }
                int intersecoes=0;
                for(int k=0;k<NumEstB;k++){
                    if(DisPtEstB[k]==Menordis){
                        SelEstB[k]=k+1;
                        intersecoes++;
                    }
                }
                if(intersecoes==1){
                    for(int k=0;k<NumEstB;k++){
                        if(SelEstB[k]>0){
                            cenario[i][j]=SelEstB[k]-1;
                        }
                    }
                }else{
                    Random gerador = new Random(intersecoes);
                    int Sorteio=(gerador.nextInt(intersecoes)+1);
                    int contador2=0;
                    int k=0,seletor=0;
                    while(contador2<Sorteio){
                        if(SelEstB[k]>0){
                            contador2++;
                            seletor=k;
                        }
                        k++;
                    }
                    cenario[i][j]=SelEstB[seletor]-1;
                }
                SelEstB=new int[NumEstB];
                
            }
        }
        }
        
        if(parmetodo==0){
            return cenario;
        }else{
            return PosEstB;
        }
        
    }
    private int[][] CenarioBase(){
        int [][] Cenaricli=new int [n*N][n*N];
        int clientemat=0,parclientemw=0;
        double number;
        for(int i=0;i<(n*N);i++){
            for(int j=0;j<(n*N);j++){
                number=(j-Math.ceil(j/n)*n)+n*(i-n*Math.ceil(i/n))+n*n*(Math.ceil(j/n)+N*Math.ceil(i/n));
                Cenaricli[i][j]=(int) number;
            }
        }
        int[]Imprimir=new int[N*n];
        for(int i=0;i<(n*N);i++){
            for(int j=0;j<(n*N);j++){
                Imprimir[j]=Cenaricli[i][j];
            }
        }
        return Cenaricli;
        
    }
    
    private int[][] Hexagono(double raio){
        double EstBdis,EstBdismeioL,aux1,aux2=l+1;
        double NumeroEsp=0;
        double cIn = 0;
        EstBdis= 2*raio*Math.cos((Math.PI/6));
        EstBdismeioL=EstBdis*Math.cos((Math.PI/6));
        if(EstBdis > n*N*l){
            EstBdis = n*N*l; 
        }
        if(EstBdismeioL > n*N*l){
            EstBdismeioL = n*N*l;
        }
        int NumeroEsp_Ml = (int) Math.round(EstBdismeioL/l);
        int NumeroEsp_ED = (int) Math.round(0.5*EstBdis/l);
        int iniFiCent = declaraInicial (2*NumeroEsp_Ml, ((N*n/2)-1));
        int iniFiVert = declaraInicial (NumeroEsp_ED, ((N*n/2)-1));        
        int pIniFiNCent = ((N*n/2)-1) + NumeroEsp_Ml;
        
        int iniFiNCent = -1;
        
        if (pIniFiNCent >= N*n){
            pIniFiNCent = ((N*n/2)-1) - NumeroEsp_Ml;
        }
        if (pIniFiNCent >= 0){
            iniFiNCent = declaraInicial (2*NumeroEsp_Ml, (pIniFiNCent));
        }
        
        int [] central = new int [N*n];
        int [] nCentral = new int [N*n];
        
        central = construirFileira (NumeroEsp_Ml,iniFiCent);
        
        if(iniFiNCent != -1) {
            nCentral = construirFileira (NumeroEsp_Ml,iniFiNCent);
        }
        System.out.println("Verifica Funcao");
        System.out.println(Arrays.toString(nCentral));
        System.out.println("Verifica Funcao2");
        System.out.println(Arrays.toString(central));
        System.out.println("Verifica Matriz");
        int [][] cobertura= new int [n*N][n*N];
        int aux_for=0;
        int aux2_for=((n*N)/2)-1;
        for (int i=iniFiVert; i<=((n*N)/2)-1;i+=NumeroEsp_ED){
            for (int j=0; j < N*n; j++){
                if(aux_for==0){
                cobertura[aux2_for][j]=central[j];
                }else{
                cobertura[aux2_for][j]=nCentral[j];
                }
            }
            aux2_for=aux2_for-NumeroEsp_ED;
            if(aux_for==0){
                aux_for=1;
            }else{
                aux_for=0;
            }
        }
        aux_for = 0;
        for (int i=(N*n/2)-1; i<N*n;i+=NumeroEsp_ED){
            for (int j=0; j < N*n; j++){
                if(aux_for==0){
                cobertura[i][j]=central[j];
                }else{
                cobertura[i][j]=nCentral[j];
                }
            }
            if (aux_for==0){
                aux_for= 1;
            }else {
                aux_for= 0;
            }
        }
        
        //Para debug
        int [] vetorDeImpressao = new int[cobertura[0].length];
        for (int i=0; i<cobertura.length;i++){
            for (int j=0; j<cobertura[0].length; j++){
            vetorDeImpressao[j] = cobertura [i][j];
        }
        System.out.println(Arrays.toString(vetorDeImpressao));
    }
        
        int teste = contaCelulas(cobertura);
        
        int [][]matrizDeCobertura = matrizDePosicoes(cobertura, teste, SR);
    /*    
        int [] vetorDeImpressao = new int[matrizDeCobertura[0].length];
        for (int i=0; i<matrizDeCobertura.length;i++){
            for (int j=0; j<matrizDeCobertura.length; j++){
            vetorDeImpressao[j] = matrizDeCobertura [i][j];
        }
        System.out.println(Arrays.toString(vetorDeImpressao));
    }
    
        
        int teste1 = contaCelulas(cobertura);
        
        System.out.println("Nova Matriz");
        int [] vetorDeImpressao1 = new int[matrizDeCobertura[0].length];
        for (int i=0; i<matrizDeCobertura.length;i++){
            for (int j=0; j<matrizDeCobertura.length; j++){
            vetorDeImpressao1[j] = matrizDeCobertura [i][j];
        }
        System.out.println(Arrays.toString(vetorDeImpressao1));
    }  
    */
        
        /*
        int NumEstB=(int) (Math.pow(Math.ceil(DisUniEst/2.0),2)+Math.pow(Math.floor(DisUniEst/2.0),2));
        //System.out.println("numero de estações base:"+NumEstB);
        int[][]cenario=new int[(n*N)][(n*N)];
        int NEBI=(int) Math.ceil(DisUniEst/2.0);
        int NEBP=(int) Math.floor(DisUniEst/2.0);
        int posini=1,contador1=1;
        
        //System.out.println(NumeroEsp);
        while(posini>=0){
            posini=(int) (((n*N)/2.0)-1-contador1*NumeroEsp);
            contador1++;
        }
        posini=(int) (posini+NumeroEsp);
        int poslin=posini;
        int poscol=posini;
        int parnp=0;
        for(int i=0;i<NumEstB;i++){
            cenario[poslin][poscol]=1;
            if((poscol+2*NumeroEsp<n*N)){
                poscol=(int) (poscol+2*NumeroEsp);
            }else{
                //i=i-1;
                if(parnp==0){
                    parnp=1;
                }else{
                    parnp=0;
                }
                poscol=(int) (posini+parnp*NumeroEsp);
                if((poslin+NumeroEsp<n*N)){
                    poslin=(int) (poslin+NumeroEsp);
                }
            }
        }
        int contador=0;
        int [][]PosEstB=new int[NumEstB][2];
        for(int i=0;i<n*N;i++){
            for(int j=0;j<n*N;j++){
                if(cenario[i][j]==1){
                    cenario[i][j]=contador;
                    PosEstB[contador][0]=i;
                    PosEstB[contador][1]=j;
                    contador++;
                }
            }
        }
        if(parmetodo==0){
        double[] DisPtEstB=new double[NumEstB];
        int[] SelEstB=new int[NumEstB];
        for(int i=0;i<n*N;i++){
            for(int j=0;j<n*N;j++){
                for(int k=0;k<NumEstB;k++){
                    DisPtEstB[k]=Math.sqrt(Math.pow((PosEstB[k][1]-j),2)+Math.pow((PosEstB[k][0]-i),2));
                }
                
                double Menordis=DisPtEstB[0];
                for(int k=0;k<NumEstB;k++){
                    if(DisPtEstB[k]<Menordis){
                        Menordis=DisPtEstB[k];
                    }
                }
                int intersecoes=0;
                for(int k=0;k<NumEstB;k++){
                    if(DisPtEstB[k]==Menordis){
                        SelEstB[k]=k+1;
                        intersecoes++;
                    }
                }
                if(intersecoes==1){
                    for(int k=0;k<NumEstB;k++){
                        if(SelEstB[k]>0){
                            cenario[i][j]=SelEstB[k]-1;
                        }
                    }
                }else{
                    Random gerador = new Random(intersecoes);
                    int Sorteio=(gerador.nextInt(intersecoes)+1);
                    int contador2=0;
                    int k=0,seletor=0;
                    while(contador2<Sorteio){
                        if(SelEstB[k]>0){
                            contador2++;
                            seletor=k;
                        }
                        k++;
                    }
                    cenario[i][j]=SelEstB[seletor]-1;
                }
                SelEstB=new int[NumEstB];
                
            }
        }
        }
        
        if(parmetodo==0){
            //return cenario;
        }else{
            //return PosEstB;
        }
        */
        return matrizDeCobertura;
    }
    
    private int declaraInicial(int disH,int pIni){
        int cIn = 0;
        int i_while1 = 0;
        while (cIn>=0){
            cIn = ((pIni))-(disH*i_while1);
            i_while1++;
        }
        cIn=cIn+(disH);
        return cIn;
    }
    private int[] construirFileira (int disH,int inifileira){
        int cIn[] = new int[n*N];
        int i_while1 = inifileira;
        while (i_while1<N*n){
            cIn[i_while1] = 1;
            i_while1=i_while1+2*disH;
        }
        return cIn;
    }
    
    private int contaCelulas (int[][] cobertura){
        int contador = 0;
        for (int i=0; i<cobertura[0].length ;i++){
            for (int j=0; j<cobertura.length; j++){
                if(cobertura[i][j]==1){
                contador+=1;
                }
            }
        }
        return contador;
    }
    
    private int [][] matrizDePosicoes (int [][] cobertura, int contador, int SR){
        int [][]matDePos = new int [contador][3];
        int aux=0;
        int hor = 0;
        int vert = 0;
        switch (SR){
            case 8:
                hor = 2;
                vert = 4;
            case 16:
                hor = 4;
                vert = 4;
            break;
            case 32:
                hor = 8;
                vert = 4;
            break;
            case 64:
                hor = 8;
                vert = 8;
            break;
            default:
        }
        for (int i=0; i<cobertura[0].length ;i++){
            for (int j=0; j<cobertura.length; j++){
                if(cobertura[i][j]==1){
                matDePos[aux][0]=i;
                matDePos[aux][1]=j;
                aux++;
                }
            }
            
        }
        aux = matDePos[0][0];
        int aux2 = 0;
        for(int i=0; i < matDePos.length; i++){
            if (matDePos[i][0]!= aux){
                aux2=i;
                break;
            }
        }
        int aux3 = 0;
        for(int i=0; i < matDePos.length; i++){
            if (matDePos[i][0]== aux){
                aux3++;
            }else{
                break;
            }
        }
        int aux4 = 0;
        for(int i=aux2; i < matDePos.length; i++){
            if (matDePos[i][0]== matDePos[aux2][0]){
                aux4++;
            }else{
                break;
            }
        }
        int aux5 = 1;
        for(int i=0; i < matDePos.length-1 ; i++){ //Generalizar o matDePos[][] para todas as colunas
            if (matDePos[i+1][0]!= matDePos[i][0]){
                aux5++;
                aux5 = aux5;
            }
        }
        System.out.println("comeca aqui");
        int [][] matSimplificada = new int [aux5][aux3+aux4];
        for (int i=0; i<aux5 ;i+=2){
            for (int j=0; j<aux3+aux4; j+=2){
                if(aux3>=aux4){
                    matSimplificada[i][j] = 1;
                    if(j+1<aux3+aux4&&i+1<aux5){
                    matSimplificada[i+1][j+1] = 1;
                    }
                }else{
                    if(j+1<aux3+aux4){
                    matSimplificada[i][j+1] = 1;
                    }
                    if(i+1<aux5){
                    matSimplificada[i+1][j] = 1;
                    }
                }
            }
        }
        
        matSimplificada = MatSimplificada(matSimplificada);
        
        int cont2=0;
        for(int i=0;i<matSimplificada.length;i++){
            for (int j=0; j<matSimplificada[0].length;j++){
                if(matSimplificada[i][j]>0){
                    matDePos[cont2][2] = matSimplificada[i][j];
                    cont2++;
                }
            }
        }
        
        int [] vetorDeImpressao2 = new int[matSimplificada[0].length];
        for (int i=0; i<matSimplificada.length;i++){
            for (int j=0; j<matSimplificada[0].length; j++){
            vetorDeImpressao2[j] = matSimplificada [i][j];
            }
        System.out.println(Arrays.toString(vetorDeImpressao2));
        }
        System.out.println("termina aqui2");
        
        return matDePos;
    }
     
    private int[][] MatSimplificada (int[][] matSimplificada){
        int delN = 16,  delH, delV, vIndex = 1, vIndexReset=1, auxH=0;
        if (SR ==16){
        delH=4;
        delV=4;
        } else if(SR==32){
        delH=4;
        delV=8;
        } else if (SR ==64){
        delH=8;
        delV=8;
        } else if (SR == 8){
        delH=2;
        delV=4;
        } else {
        delH=0;
        delV=0;
        }
        for (int i = 0; i<matSimplificada.length;i++){
            for (int j = 0; j<matSimplificada[0].length; j++){
                if (matSimplificada[i][j] == 1){
                    matSimplificada[i][j] = vIndex;
                    auxH++;
                }
                if (auxH == delH){
                    auxH=0;
                    vIndex++;
                }
            }
            auxH=0;
            if((i+1)%delV==0){
                vIndex++;
                vIndexReset = vIndex;
            } else {
                vIndex = vIndexReset;
            }
        }
        return matSimplificada;
    }
}


