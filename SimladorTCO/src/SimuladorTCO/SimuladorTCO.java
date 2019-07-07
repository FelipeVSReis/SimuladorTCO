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
public class SimuladorTCO{
    //parametros do modelo manhatttan
    int N=4; //Número de Blocos numa fileira
    int n=4; //Número de Resisdencias numa fileira
    double l=0.03; // Dsitancia entre residencias (km)
    double L=n*l;        // Distancia Entre Blocos
    double num_andares=1;
    //parametros de custos
    double sal=190; //($/h)
    double SLAR=10; // Multa Residencial($/h)
    double SLAC=100; // Multa Residencial($/h)
    double PrSpleaYe=204;//($/Ye) ALuguel de Espectro por Link MW
    double alpha=1.1; //fator de impacto
    double PrKw=0.1; // Preço Energia
    double PrGbsh=3; // Preço Aluguel Taxa MW
    double paraluext=1.3;//parametro eevaror de aluiguel MW
    //double FatorMultiplicador=10;
    double vel=20; //(Km/h)
    double PorcentagemCom=0.8; // porcentagem usuarios comercias
    int SR=32;   //splitting ratio(16,32 ou 64)
    int NPC=72;  //numero de portas olt por chassi botar 72
   //subtituir e arqpad esqprot
   //seletores
    int selCO=1;
    // 1 olt 
    // 2 olt cruzado
    int selENA=1;
    // 1 fibra
    // 2 fibra duplicada
    int selRN=1;
    // 1 splitter
    // 2 awg
    int selEND=1;
    // 1 fibra
    // 2 fibra duplicada
    int selRES=1;
    // 1 ONU
    // 2 ONU5g
    int TypeMW=0;
    // 0 sem MW
    // 1 MW Leas.
    // 2 MW upgrade
    // 3 MW Prop.
    double PSFP=1; 
    double P_low=37;
    double P_high=95.2;

    //Parametros Monte Carlo
    //int numittt=200;
    int Ny=5;
    int Num_Tentativas=1000;
    double PorcRep=0.3; //reparo
    //Dados de Equipamentos
    // {onu,fiber,splitter,olt port,Rn chassi,Olt Chassi,switches,GES,Micro,Pico,Antena,macro}
    double[] taxas_de_falha={256,2381,120,1075,666.6,500,200,500,3225.8,1612.9,540,32258.06}; //(fit)
    double[] tempo_de_reparo={1,7,1,1,1,1,2,1,2,2,1,7}; //(h)
    double[] tempo_de_inst={1,0,0.166666667,0.166666667,0.166666667,0.5,0.166666667,1,2,1,0.16,24}; //(h)
    double[] preco={350,0,50,7600,700,4500,367*Math.ceil((num_andares/16.0)),1980,9000,1600,2000,22000}; //($)
    double[] KWh={5,0,0,1197,0,0,0,20+14*num_andares,150,45,0,22000}; //(W)
    double[] areacobertura={0.25,0.1};//{micro,pico}//até 0.25
    double PkmVF=160; //instalar
    double PkmTF=130000; //trenching
    
    
    public SimuladorTCO()// throws IOException
    { 
    if(TypeMW==0){
        taxas_de_falha[8]=0;
        taxas_de_falha[9]=0;
        taxas_de_falha[10]=0;
        taxas_de_falha[11]=0;
        KWh[8]=0;
        KWh[9]=0;
        KWh[10]=0;
        KWh[11]=0;
    }
    //System.out.println(Arrays.toString(taxas_de_falha));
    //Numero de Equipamentos
    //{onu,dstep,spt,rn,fstep,oltc,fswt,pfstep,dswt,pdstep,GES,Micro,Pico,Ant,oltp,macro, olt splitter, olt splitter switch}
    
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
    //SELETORES
    int suboltspl=0;
    if(selCO>0){
        numspq=(int) Math.ceil((double) Math.pow(n,2)/SR);
        numEq[14]=((int) (Math.pow(N,2)*numspq));
        numEq[5]=(int) Math.ceil((double) ((int) (Math.pow(N,2)*numspq))/NPC);
        if(selCO==2){
            if(numEq[5]%2==0){
                suboltspl=(int) ((Math.ceil(((numEq[14]*1.0)/NPC))*NPC)-numEq[14]);
                //System.out.println("entrou no 1" + suboltspl);
            }
            else if(numEq[5]%2!=0){
                suboltspl=numEq[14]-((numEq[5]-1)*NPC);
            }
            numEq[16]=numEq[14]-suboltspl;
            numEq[17]=numEq[14];
        }
    }
    if(selENA>0){
    numEq[4]=(int) ((N*(N+1)));
        if(selENA==2){
            numspq=(int) Math.ceil((double) Math.pow(n,2)/SR);
            numEq[6]=((int) (Math.pow(N,2)*numspq));
            numEq[7]=(N-1)*N;
        }
    }
    if(selRN>0){
        numEq[3]=(int) Math.pow(N,2);
        numspq=(int) Math.ceil((double) Math.pow(n,2)/SR);
        numEq[2]=((int) (Math.pow(N,2)*numspq));
        if(selEND==2){ 
            //por enquanto o mesmo
            numEq[3]=(int) Math.pow(N,2);
            numspq=(int) Math.ceil((double) Math.pow(n,2)/SR);
            numEq[2]=((int) (Math.pow(N,2)*numspq));
        }
    }
    if(selEND>0){
    numEq[1]=(int) ((n*(n+1))*Math.pow(N,2));
        if(selEND==2){
            numEq[8]=n*n*N*N;
            numEq[9]=(n-1)*n*N*N;
        }
    }
    if(selRES>0){
        numEq[0]=(int) (Math.pow(n,2)*Math.pow(N,2));
        if(selRES==2){
            //Por Equato o mesmo
            numEq[0]=(int) (Math.pow(n,2)*Math.pow(N,2));
        }
    }
    CenarioMWCli=new int[2][2];
    CenarioMWPEB=new int[2][2];
    if(TypeMW==3){
        CenarioMWCli=MatrizEstBase(areacobertura[1], 0);
        CenarioMWPEB=MatrizEstBase(areacobertura[1], 1);
        numEq[12]=CenarioMWPEB.length;
        numEq[13]=(int) Math.ceil((numEq[12]/16.0));
        numEq[15]=1;
    }
    
    int numtoteq=0;
    for(int i=0;i<18;i++){
        numtoteq=numtoteq+numEq[i];
    }
    int NumeroOLTPextra=(int) ((Math.ceil(((numEq[2]*1.0)/NPC))*NPC)-numEq[2]);
    int NumeroOLTPextra2=numEq[2]-((numEq[5]-1)*NPC);
    if(((numEq[5]%2==0)&&((numEq[2]/NPC) %2!=0))&&selCO==2){
        numtoteq=numtoteq+NumeroOLTPextra;
    }
    if((numEq[5]%2!=0)&&selCO==2){
        numtoteq=numtoteq+1+NumeroOLTPextra2;
    }
    //System.out.println("{onu,dstep,spt,rn,fstep,oltc,fswt,pfstep,dswt,pdstep,GES,Micro,Pico,Ant,oltp,macro, olt splitter, olt splitter switch}");
    //System.out.println(Arrays.toString(numEq));
    //System.out.println("numeroeq" + numtoteq);
    // Modelo Geometrico 
    //criando vetor eqipamentos
    
    
    Equipamento [] Equi = new Equipamento[numtoteq];
    for(int i=0;i<numtoteq;i++){
        Equi[i] = new Equipamento();
    }
    int parpas=0;
    // Arquitetura Principal
    //Num Eq
    /*
    Onu-0
    FDS-1
    Splitters-2
    Olt ports-3
    Rn_Chassi-4
    FFS-5
    Olt_Chassi-6
    Olt_Switch-7
    PFFS-8
    Olt_Switch-9
    PDFS-10
    GES-11
    Macro-12
    Antena_Macro-13
    Micro-14
    Pico-15
    */
    //numero médio de falhas
    int NumEqNP=numEq[0]+numEq[1]+numEq[2]+numEq[2]+numEq[3]+numEq[4]+numEq[5];
    //System.out.println(NumEqNP);
    int Nummedclifal=(int) Math.ceil(((N*N*n*n*(14+n+N))/(2.0*NumEqNP)));
    System.out.println("med"+Nummedclifal);
    
    //Onus     
    int [][] Onucli=new int[3][numEq[0]]; //cria armazenador de posição de ONU
    if(numEq[0]>0){
        int paronu=0;
        //preco,Qaux,PVaux,PHaux,Cli,Cli,Dis,TXf,TemR,TemI.
        for(int i=0;i<numEq[0];i++){
            double[] ONU=CriarONU(i,paronu);
            Equi[parpas].custoeq=ONU[0];
            Onucli[0][i]=(int) ONU[1];
            Onucli[1][i]=(int) ONU[2];
            Onucli[2][i]=(int) ONU[3];
            Equi[parpas].clientes[0]=(int) ONU[4];
            Equi[parpas].clientes[1]=(int) ONU[5];
            Equi[parpas].distancia=ONU[6];
            Equi[parpas].taxadefalha=ONU[7];
            Equi[parpas].tempodereparo=ONU[8];
            Equi[parpas].instalacao=ONU[9];
            Equi[parpas].ConsumoEnergia=KWh[0];
            Equi[parpas].tipo=0;
            paronu++;
            if(paronu==Math.pow(n,2)){
                paronu=0;
            }
            parpas++;
        }
    }
    //System.out.println(parpas);
    //Fiber distribuition Steps
    if(numEq[1]>0){
        int pards=0;
        for(int i=0;i<numEq[1];i++){
            double[] DS=CriarDisStep(i,parpas,pards,numEq[0],Onucli);
            //Cli,Cli,Dis,Prot,Prot,Prot,TXf,TemR.
            Equi[parpas].clientes[0]=(int) DS[0];
            Equi[parpas].clientes[1]=(int) DS[1];
            Equi[parpas].distancia=DS[2];
            Equi[parpas].protecao[0]=(int) DS[3];
            Equi[parpas].protecao[1]=(int) DS[4];
            Equi[parpas].protecao[2]=(int) DS[5];
            //System.out.println(Arrays.toString(DS));
            Equi[parpas].taxadefalha=DS[6];
            Equi[parpas].tempodereparo=DS[7];
            Equi[parpas].tipo=1;
            pards++;
            if(pards==((n+1)*n)){
                pards=0;
            }
            parpas++;
        }
    }
    //Splitters
    if(numEq[2]>0){
        int Cps=0;
        int CliQ=(int) Math.pow(n,2);
        for(int i=0;i<numEq[2];i++){
            double [] Splitter=Criar_Splitter(i, parpas, numspq, numEq[0], Cps, CliQ);
            //Cli,Cli,Dis,preco,TXf,TemR,TemI,Cps.
            Equi[parpas].clientes[0]=(int) Splitter[0];
            Equi[parpas].clientes[1]=(int) Splitter[1];
            Equi[parpas].distancia=Splitter[2];
            Equi[parpas].custoeq=Splitter[3];
            Equi[parpas].taxadefalha=Splitter[4];
            Equi[parpas].tempodereparo=Splitter[5];
            Equi[parpas].instalacao=Splitter[6];
            Equi[parpas].instalacao=Splitter[6];
            Equi[parpas].tipo=2;
            Cps=(int) Splitter[7];
            CliQ=(int) Splitter[8];
            parpas++;
        }
    }
    //OLT ports
    if(numEq[2]>0){
        int Cps=0;
        int CliQ=(int) Math.pow(n,2);
        for(int i=0;i<numEq[2];i++){
            double [] Splitter=Criar_Splitter(i, parpas, numspq, numEq[0], Cps, CliQ);
            //Cli,Cli,Dis,preco,TXf,TemR,TemI,Cps.
            Equi[parpas].clientes[0]=(int) Splitter[0];
            Equi[parpas].clientes[1]=(int) Splitter[1];
            Equi[parpas].taxadefalha=taxas_de_falha[3]/1000000000.0;
            Equi[parpas].tempodereparo=tempo_de_reparo[3];
            Equi[parpas].instalacao=tempo_de_inst[3];
            Equi[parpas].custoeq=preco[3];
            Equi[parpas].ConsumoEnergia=KWh[3];
            Equi[parpas].tipo=3;
            Cps=(int) Splitter[7];
            CliQ=(int) Splitter[8];
            parpas++;
        }
        if(((numEq[5]%2==0)&&((numEq[2]/NPC) %2!=0))&&selCO==2){
            for(int i=0;i<NumeroOLTPextra;i++){
                Equi[parpas].clientes[0]=0;
                Equi[parpas].clientes[1]=0;
                Equi[parpas].taxadefalha=taxas_de_falha[3]/1000000000.0;
                Equi[parpas].tempodereparo=tempo_de_reparo[3];
                Equi[parpas].instalacao=tempo_de_inst[3];
                Equi[parpas].custoeq=preco[3];
                Equi[parpas].ConsumoEnergia=KWh[3];
                Equi[parpas].tipo=3;
            parpas++;
            }
        }else if((numEq[5]%2!=0)&&selCO==2){
            for(int i=0;i<NumeroOLTPextra2;i++){
                Equi[parpas].clientes[0]=0;
                Equi[parpas].clientes[1]=0;
                Equi[parpas].taxadefalha=taxas_de_falha[3]/1000000000.0;
                Equi[parpas].tempodereparo=tempo_de_reparo[3];
                Equi[parpas].instalacao=tempo_de_inst[3];
                Equi[parpas].custoeq=preco[3];
                Equi[parpas].ConsumoEnergia=KWh[3];
                Equi[parpas].tipo=3;
            parpas++;
            }
        }
        if(selCO==2){
            int extra;
            if(numEq[5]%2!=0){
                extra=NumeroOLTPextra2;
            }else{
                extra=NumeroOLTPextra;
            }
            int NCtrab;
            int Nolttot=numEq[14]+extra;
            int [] Vetordef=new int[Nolttot];
            if(numEq[5]%2==0){
                NCtrab=numEq[5];
            }else{
                NCtrab=numEq[5]+1;
            }
            int varsai=0;
            int countwhile=0;
            int countwhile2=0;
            while(varsai!=1){
                if(Equi[countwhile].tipo==3){
                    Vetordef[countwhile2]=countwhile;
                    countwhile2++;
                    if(countwhile2==Nolttot){
                        varsai=1;
                    }
                }
                countwhile++;
            }
            int [] confoltchassi=new int[NCtrab];
            //int confoltchassiprot=0;
            int dstoltp=numEq[14];
            for(int i=0;i<NCtrab;i++){
                if(numEq[5]%2!=0){
                    if((dstoltp-NPC)<0){
                        confoltchassi[i]=dstoltp;
                        confoltchassi[i+1]=dstoltp;
                        i++;
                    }else{
                        confoltchassi[i]=NPC;
                        dstoltp=dstoltp-NPC;
                    }
                }else{
                    confoltchassi[i]=NPC; 
                }
            }
            //comecar daqui
            for(int i=1;i<NCtrab;i++){
                confoltchassi[i]=confoltchassi[i]+confoltchassi[i-1];
            }
            int seletoresp=0;
            for(int i=0;i<Nolttot;i++){
                if(seletoresp%2==0){
                    if(i<confoltchassi[seletoresp]){
                        int dif;
                        if(seletoresp==0){
                            dif=confoltchassi[0];
                        }else{
                            dif=confoltchassi[seletoresp]-confoltchassi[seletoresp-1];;
                        }
                        Equi[Vetordef[i]].protecao[0]=8;
                        Equi[Vetordef[i]].protecao[1]=Equi[Vetordef[i]+dif].clientes[0];
                        Equi[Vetordef[i]].protecao[2]=Equi[Vetordef[i]+dif].clientes[1];
                       }else{
                        seletoresp++;
                        int dif=confoltchassi[seletoresp]-confoltchassi[seletoresp-1];
                        Equi[Vetordef[i]].protecao[0]=8;
                        Equi[Vetordef[i]].protecao[1]=Equi[Vetordef[i]-dif].clientes[0];
                        Equi[Vetordef[i]].protecao[2]=Equi[Vetordef[i]-dif].clientes[1];
                    }
                }else{
                    if(i<confoltchassi[seletoresp]){
                        int dif=confoltchassi[seletoresp]-confoltchassi[seletoresp-1];;
                        Equi[Vetordef[i]].protecao[0]=8;
                        Equi[Vetordef[i]].protecao[1]=Equi[Vetordef[i]-dif].clientes[0];
                        Equi[Vetordef[i]].protecao[2]=Equi[Vetordef[i]-dif].clientes[1];
                       }else{
                        seletoresp++;
                        int dif=confoltchassi[seletoresp]-confoltchassi[seletoresp-1];
                        Equi[Vetordef[i]].protecao[0]=8;
                        Equi[Vetordef[i]].protecao[1]=Equi[Vetordef[i]+dif].clientes[0];
                        Equi[Vetordef[i]].protecao[2]=Equi[Vetordef[i]+dif].clientes[1];
                    }
                }
            }
            for(int i=0;i<Nolttot;i++){
                if(Equi[Vetordef[i]].protecao[1]== Equi[Vetordef[i]].protecao[2]){
                    Equi[Vetordef[i]].protecao[0]=9;//change value
                }
            }
        }
    }
    //Rn Chassi
    if(numEq[3]>0){
        int CpRn=0;
        int RCliQ=(int) Math.pow(n,2);
        for(int i=0;i<numEq[3];i++){
            double [] Rn_Chassi=CriarRn_Chassi(i, CpRn, RCliQ, numEq[0]);
            //Cli,Cli,Dis,preco,TXf,TemR,TemI,CpRn.
            Equi[parpas].clientes[0]=(int) Rn_Chassi[0];
            Equi[parpas].clientes[1]=(int) Rn_Chassi[1];
            Equi[parpas].distancia=Rn_Chassi[2];
            Equi[parpas].custoeq=Rn_Chassi[3];
            Equi[parpas].taxadefalha=Rn_Chassi[4];
            Equi[parpas].tempodereparo=Rn_Chassi[5];
            Equi[parpas].instalacao=Rn_Chassi[6];
            Equi[parpas].tipo=4;
            CpRn=(int) Rn_Chassi[7];
            parpas++;
        }
    }
    //Quadras
    int NQ=(int) Math.pow(N,2);
    int [][] Quadracli = new int[2][NQ];
    for(int i=0;i<NQ;i++){
        int PosAuxFFV=(int) Math.floor((double) i/N);
        int PosAuxFFH= i-N*PosAuxFFV;
        Quadracli[0][i]=PosAuxFFV;
        Quadracli[1][i]=PosAuxFFH;
    }
    //Fiber feder Steps
    if(numEq[4]>0){
        for(int i=0;i<numEq[4];i++){
            double[] FS=CriarFedStep(i, parpas, Quadracli, numEq[0], NQ);
            //System.out.println(Arrays.toString(FS));
            //Cli,Cli,Dis,Prot,Prot,Prot,TXf,TemR.
            Equi[parpas].clientes[0]=(int) FS[0];
            Equi[parpas].clientes[1]=(int) FS[1];
            Equi[parpas].distancia=FS[2];
            Equi[parpas].protecao[0]=(int) FS[3];
            Equi[parpas].protecao[1]=(int) FS[4];
            Equi[parpas].protecao[2]=(int) FS[5];
            Equi[parpas].taxadefalha=FS[6];
            Equi[parpas].tempodereparo=FS[7];
            Equi[parpas].tipo=5;
            parpas++;
        }
    }
    //olt chassi
    if(numEq[5]>0){
        int alnumspt=numEq[2];
        int contoc = 0;
        for(int i=0;i<numEq[5];i++){
            double[] OLTC=CriarOltC(i, alnumspt, numEq[0], numEq[1], numEq[2], contoc, Equi);
            //Cli,Cli,Preco,TXf,TemR,TemI,contoc.
            Equi[parpas].clientes[0]=(int) OLTC[0];
            Equi[parpas].clientes[1]=(int) OLTC[1];
            Equi[parpas].custoeq=OLTC[2];
            Equi[parpas].taxadefalha=OLTC[3];
            Equi[parpas].tempodereparo=OLTC[4];
            Equi[parpas].instalacao=OLTC[5];
            Equi[parpas].tipo=6;
            contoc=(int) OLTC[6];
            parpas++;
        }
        
        if(selCO==2){
                int numtotchassi=numEq[5];
                if(numEq[5]%2!=0){
                    Equi[parpas].clientes[0]=0;
                    Equi[parpas].clientes[1]=0;
                    Equi[parpas].custoeq=Equi[parpas-1].custoeq;
                    Equi[parpas].taxadefalha=Equi[parpas-1].taxadefalha;
                    Equi[parpas].tempodereparo=Equi[parpas-1].tempodereparo;
                    Equi[parpas].instalacao=Equi[parpas-1].instalacao;
                    Equi[parpas].tipo=6;
                    parpas++;
                    numtotchassi=numtotchassi+1;
                }
                for(double i=0.5;i<numtotchassi;i++){
                    int cima=(int) Math.ceil((i*1.0));
                    int baixo=(int) Math.floor((i*1.0));
                    Equi[parpas-(numtotchassi)+cima].protecao[0]=8;
                    Equi[parpas-(numtotchassi)+cima].protecao[1]=Equi[parpas-(numtotchassi)+baixo].clientes[0];
                    Equi[parpas-(numtotchassi)+cima].protecao[2]=Equi[parpas-(numtotchassi)+baixo].clientes[1];
                    Equi[parpas-(numtotchassi)+baixo].protecao[0]=8;
                    Equi[parpas-(numtotchassi)+baixo].protecao[1]=Equi[parpas-(numtotchassi)+cima].clientes[0];
                    Equi[parpas-(numtotchassi)+baixo].protecao[2]=Equi[parpas-(numtotchassi)+cima].clientes[1];
                    i++;
                }
        }
        
        
    }
    //{onu,dstep,spt,rn,fstep,oltc,fswt,pfstep,dswt,pdstep,GES,Micro,Pico,Ant,oltp,macro}
    //olt swithchs
    if(numEq[6]>0){
        int Cps=0;
        int CliQ=(int) Math.pow(n,2);
        for(int i=0;i<numEq[2];i++){
            double [] Splitter=Criar_Splitter(i, parpas, numspq, numEq[0], Cps, CliQ);
            //Cli,Cli,Dis,preco,TXf,TemR,TemI,Cps.
            Equi[parpas].clientes[0]=(int) Splitter[0];
            Equi[parpas].clientes[1]=(int) Splitter[1];
            Equi[parpas].taxadefalha=taxas_de_falha[6]/1000000000.0;
            Equi[parpas].tempodereparo=tempo_de_reparo[6];
            Equi[parpas].instalacao=tempo_de_inst[6];
            Equi[parpas].custoeq=preco[6];
            Equi[parpas].tipo=7;
            Cps=(int) Splitter[7];
            CliQ=(int) Splitter[8];
            parpas++;
        }
    }
    //feeder fiber protection steps
    if(numEq[7]>0){
        for(int i=0;i<numEq[7];i++){
            double[] PFFS=CriarPFedStep(i);
            //prot,prot,prot,dis,TXf,TemR.
            Equi[parpas].protecao[0]=(int) PFFS[0];
            Equi[parpas].protecao[1]=(int) PFFS[1];
            Equi[parpas].protecao[2]=(int) PFFS[2];
            Equi[parpas].distancia= PFFS[3];     
            Equi[parpas].taxadefalha=PFFS[4];
            Equi[parpas].tempodereparo=PFFS[5];
            Equi[parpas].tipo=8;
            parpas++;
        }
    }
    //onu switch
    if(numEq[8]>0){
        for(int i=0;i<numEq[0];i++){
            Equi[parpas].clientes=Equi[i].clientes;
            //Equi[parpas].distancia=Equi[i].distancia;
            int Quadraaux=(int) Math.floor((double) i/Math.pow(n,2));
            int PosAuxFFV=(int) Math.floor((double) Quadraaux/N);
            int PosAuxFFH= Quadraaux-N*PosAuxFFV;
            if(PosAuxFFV<N/2.0){
                if(PosAuxFFH<N/2.0){
                    Equi[parpas].distancia=(((N-(2*PosAuxFFH+1))+(N-(2*PosAuxFFV+1)))*(L/2.0));
                }else{
                    Equi[parpas].distancia=(((N-(2*(N-PosAuxFFH-1)+1))+(N-(2*PosAuxFFV+1)))*(L/2.0));
                }
            }else{
                PosAuxFFV=N-PosAuxFFV-1;
                if(PosAuxFFH<N/2.0){
                    Equi[parpas].distancia=(((N-(2*PosAuxFFH+1))+(N-(2*PosAuxFFV+1)))*(L/2.0));
                }else{
                    Equi[parpas].distancia=(((N-(2*(N-PosAuxFFH-1)+1))+(N-(2*PosAuxFFV+1)))*(L/2.0));
                }
            }
            Equi[parpas].taxadefalha=taxas_de_falha[6]/1000000000.0;
            Equi[parpas].tempodereparo=tempo_de_reparo[6];
            Equi[parpas].instalacao=tempo_de_inst[6];
            Equi[parpas].custoeq=preco[6];
            Equi[parpas].tipo=9;
            parpas++;
        }
    }
    //distribuition fiber protection steps
    if(numEq[9]>0){
        int pardsp=0;
        for(int i=0;i<numEq[9];i++){
            double[] PFDS=CriarPdisStep(i, pardsp);
            //prot,prot,prot,dis,TXf,TemR.
            Equi[parpas].protecao[0]=(int) PFDS[0];
            Equi[parpas].protecao[1]=(int) PFDS[1];
            Equi[parpas].protecao[2]=(int) PFDS[2];
            Equi[parpas].distancia=PFDS[3];
            Equi[parpas].taxadefalha=PFDS[4];
            Equi[parpas].tempodereparo=PFDS[5];
            Equi[parpas].tipo=10;
            pardsp++;
            if(pardsp==((n-1)*n)){
                pardsp=0;
            }
            parpas++;
        }
    }
    //GES
    if(numEq[10]>0){
        int paronu=0;
        for(int i=0;i<numEq[10];i++){
            double[] ONU=CriarONU(i,paronu);
            Equi[parpas].clientes[0]=i;
            Equi[parpas].clientes[1]=i;
            Equi[parpas].distancia=ONU[6];
            Equi[parpas].taxadefalha=taxas_de_falha[7]/1000000000.0;
            Equi[parpas].tempodereparo=tempo_de_reparo[7];
            Equi[parpas].instalacao=tempo_de_inst[7];
            Equi[parpas].custoeq=preco[7];
            Equi[parpas].ConsumoEnergia=KWh[7];
            Equi[parpas].tipo=11;
            paronu++;
            if(paronu==Math.pow(n,2)){
                paronu=0;
            }
            parpas++;
        }
        
    }
    //Macro
    if(numEq[15]>0){
        Equi[parpas].protecao[0]=7;
        Equi[parpas].taxadefalha=taxas_de_falha[11]/1000000000.0;
        Equi[parpas].tempodereparo=tempo_de_reparo[11];
        Equi[parpas].instalacao=tempo_de_inst[11];
        Equi[parpas].custoeq=preco[11];
        Equi[parpas].ConsumoEnergia=KWh[11];
        Equi[parpas].tipo=12;
        parpas++;
        
    }
    //Antena Macro
    if(numEq[13]>0){
        int NuEB=0;
        int contEB=0;
        if(numEq[11]>0){
            NuEB=numEq[11];
        }else if(numEq[12]>0){
            NuEB=numEq[12];
        }
        for(int i=0;i<numEq[13];i++){
        if(TypeMW>0){
            Equi[parpas].protecao[0]=7;
            Equi[parpas].protecao[1]=contEB;
            if(NuEB<=16){
                Equi[parpas].protecao[2]=contEB+NuEB-1;
            }else{
                Equi[parpas].protecao[2]=contEB+15;
                NuEB=NuEB-16;
                contEB=contEB+16;
            }
        }
        Equi[parpas].taxadefalha=taxas_de_falha[10]/1000000000.0;
        Equi[parpas].tempodereparo=tempo_de_reparo[10];
        Equi[parpas].instalacao=tempo_de_inst[10];
        Equi[parpas].custoeq=preco[10];
        Equi[parpas].ConsumoEnergia=KWh[10];
        Equi[parpas].tipo=13;
        parpas++;
        }
    }
    //Micro
    if(numEq[11]>0){
        for(int i=0;i<numEq[11];i++){
            Equi[parpas].distancia=1;
            Equi[parpas].protecao[0]=7;
            Equi[parpas].protecao[1]=i;
            Equi[parpas].protecao[2]=i;
            Equi[parpas].taxadefalha=taxas_de_falha[8]/1000000000.0;
            Equi[parpas].tempodereparo=tempo_de_reparo[8];
            Equi[parpas].instalacao=tempo_de_inst[8];
            Equi[parpas].custoeq=preco[8];
            Equi[parpas].ConsumoEnergia=KWh[8];
            Equi[parpas].tipo=14;
            parpas++;
        }
    }
    //Pico
    if(numEq[12]>0){
        for(int i=0;i<numEq[12];i++){
            Equi[parpas].distancia=1;
            Equi[parpas].protecao[0]=7;
            Equi[parpas].protecao[1]=i;
            Equi[parpas].protecao[2]=i;
            Equi[parpas].taxadefalha=taxas_de_falha[9]/1000000000.0;
            Equi[parpas].tempodereparo=tempo_de_reparo[9];
            Equi[parpas].instalacao=tempo_de_inst[9];
            Equi[parpas].custoeq=preco[9];
            Equi[parpas].ConsumoEnergia=KWh[9];
            Equi[parpas].tipo=15;
            parpas++;
        }
    }
    int numeqant=0;
    for(int i=6;i<16;i++){
        numeqant=numeqant+numEq[i];
    }
    //Olt port Switch
    if((numEq[17]>0)&&(selCO==2)){
        int Cps=0;
        int CliQ=(int) Math.pow(n,2);
        for(int i=0;i<numEq[17];i++){
            double [] Splitter=Criar_Splitter(i, parpas, numspq, numEq[0], Cps, CliQ);
            Equi[parpas].distancia=1;
            Equi[parpas].clientes[0]=(int) Splitter[0];
            Equi[parpas].clientes[1]=(int) Splitter[1];
            Equi[parpas].taxadefalha=taxas_de_falha[6]/1000000000.0;
            Equi[parpas].tempodereparo=tempo_de_reparo[6];
            Equi[parpas].instalacao=tempo_de_inst[6];
            Equi[parpas].custoeq=preco[6];
            Equi[parpas].ConsumoEnergia=KWh[6];
            Cps=(int) Splitter[7];
            CliQ=(int) Splitter[8];
            Equi[parpas].tipo=16;
            parpas++;
        }
    }
    //Olt Spliter
    
    if((numEq[16]>0)&&(selCO==2)){
        int Cps=0;
        int CliQ=(int) Math.pow(n,2);
        int extra;
        if(numEq[5]%2!=0){
            extra=NumeroOLTPextra2;
        }else{
            extra=NumeroOLTPextra;
        }
        int Nolttot=numEq[14]+extra; 
        int [] Vetordef=new int[Nolttot];
        int varsai=0;
        int countwhile=0;
        int countwhile2=0;
        while(varsai!=1){
            if(Equi[countwhile].tipo==3){
                Vetordef[countwhile2]=countwhile;
                countwhile2++;
                if(countwhile2==Nolttot){
                    varsai=1;
                }
            }
                countwhile++;
        }
        for(int i=0;i<Nolttot;i++){
            if(Equi[Vetordef[i]].protecao[0]==8){
                if(Equi[Vetordef[i]].clientes[0]==0&&Equi[Vetordef[i]].clientes[1]==0){
                }else{
                    Equi[parpas].tipo=17;
                    Equi[parpas].distancia=1;
                    Equi[parpas].clientes[0]=Equi[Vetordef[i]].clientes[0];
                    Equi[parpas].clientes[1]=Equi[Vetordef[i]].clientes[1];
                    Equi[parpas].protecao[0]=Equi[Vetordef[i]].protecao[0];
                    Equi[parpas].protecao[1]=Equi[Vetordef[i]].protecao[1];
                    Equi[parpas].protecao[2]=Equi[Vetordef[i]].protecao[2];
                    Equi[parpas].taxadefalha=taxas_de_falha[2]/1000000000.0;
                    Equi[parpas].tempodereparo=tempo_de_reparo[2];
                    Equi[parpas].instalacao=tempo_de_inst[2];
                    Equi[parpas].custoeq=preco[2];
                    Equi[parpas].ConsumoEnergia=KWh[2];
                    parpas++;
                }
            }
        }
    }
    
    /*
    System.out.println("nEquipamento:  " +numtoteq);
    for(int i=0;i<numtoteq;i++){
    System.out.println("------------------------" + (i+1));    
    System.out.println("Equipamento:  " + (i+1));
    System.out.println("Tipo:  " + Equi[i].tipo);
    System.out.println("Clientes:  "+ Arrays.toString(Equi[i].clientes));
    System.out.println("Proteção:  "+Arrays.toString(Equi[i].protecao));
    System.out.println("Distancia:  "+Equi[i].distancia  );
    System.out.println("Preço:  "+Equi[i].custoeq);
    System.out.println("Tempo de Instalação:  "+Equi[i].instalacao);
    System.out.println("Taxa de Falha:  "+Equi[i].taxadefalha);
    System.out.println("Tempo de Reparo:  "+Equi[i]. tempodereparo);
    }*/
    
    double sumaux3=0;
        double sumaux4=0;
        for(int i=0;i<(n/2.0);i++){
            for(int j=0;j<(n/2.0);j++){
                sumaux3=sumaux3+(((n-2*i-1)/2)+((n-2*j-1)/2));
            }
        }
        for(int i=0;i<(N/2.0);i++){
            for(int j=0;j<(N/2.0);j++){
                sumaux4=sumaux4+(((N-2*i-1)/2)+((N-2*j-1)/2));
            }
        }
        double Custotfiber2=(4*sumaux3*l*N*N+4*sumaux4*L)*PkmVF;
        //System.out.println(Custotfiber2);
        
    //Taxa de Falha Geral
        double TXFG=0;
        for(int i=0;i<numtoteq;i++){
            TXFG=TXFG+Equi[i].taxadefalha;
        }
        System.out.println("taxa de falha geral:  "+TXFG);
        //System.out.println(TXFG);
    //tempo e numero de iterações
        /*
        double TEF;
        TEF=1000000000/(1000000000*(TXFG+TXFG*TXFG));
        TEF=(365*24)/TEF;
        TEF=Math.ceil(TEF)*Ny*2;
        int numittt=(int) TEF;
        System.out.println("numero de iterações" + numittt);*/
    //numero de falhas de ONU , FDS e Olt Chassi em Ny
        double TEFOnu=numEq[0]*((365*24*Ny)/(1000000000.0/taxas_de_falha[0]));
        //System.out.println("Onu"+TEFOnu);
    //Custo
        double CIPG=0;
        for(int i=0;i<numtoteq;i++){
            CIPG=CIPG+Equi[i].taxadefalha*Equi[i].custoeq*PorcRep;
        }
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
                    for(int j=0;j<Num_i;j++){
                    Matriz_Parametro[0][j]=Equifal[j].clientes[0];
                    Matriz_Parametro[1][j]=Equifal[j].clientes[1];
                    Matriz_Parametro[2][j]=Equifal[j].tempodereparo;      
                    Matriz_Parametro[3][j]=Equifal[j].distancia;
                    Matriz_Parametro[4][j]=Equifal[j].protecao[0];
                    Matriz_Parametro[5][j]=Equifal[j].protecao[1];
                    Matriz_Parametro[6][j]=Equifal[j].protecao[2];
                    Matriz_Parametro[7][j]=1.0;
                    Matriz_Parametro[8][j]=Equifal[j].tipo;
                    }
                    double[] VetorcUSTO=new double[Num_i];
                    for(int j=0;j<Num_i;j++){
                    Matriz_Parametro[7][j]=0.0;    
                    VetorcUSTO[j]=ObterCusto(numEq[0],Num_i,Matriz_Parametro,PorcentagemCom, CenarioMWCli, Cenaricli);
                    
                    Matriz_Parametro[7][j]=1.0;
                    }
                    //System.out.println(Arrays.toString(VetorcUSTO));
                    
                    double CustoSLAini=VetorcUSTO[0];
                    possel=0;
                    double sumtxfds=0;
                    for(int j=0;j<Num_i;j++){
                        sumtxfds=sumtxfds+Equi[posfalha[j]].taxadefalha;
                        if(CustoSLAini<VetorcUSTO[j]){
                        CustoSLAini=VetorcUSTO[j];
                        possel=j;
                        }
                    }
                    
                    Lambda=TXFG+(1/(Equi[posfalha[possel]].tempodereparo+(Equi[posfalha[possel]].distancia/vel)))-sumtxfds;
                    double Tempo=(1/Lambda);//*(Math.log(u));
                    //não generico
                    double CustoSLA;
                    double CustoET;
                    CustoSLA=ObterCusto(numEq[0],Num_i,Matriz_Parametro,PorcentagemCom, CenarioMWCli, Cenaricli);
                    if(TypeMW==1){
                    CustoAlugueletapa=ObterAluguel(numEq[0], Num_i, Matriz_Parametro, PorcentagemCom, CenarioMWCli, Cenaricli, Nummedclifal);
                    }
                    int numfalha;
                    
                    
                    double Custoret=0;
                    /*
                    for(int j=0;j<Num_i;j++){
                        Custoret=Custoret+((Equi[posfalha[j]].taxadefalha+Equi[posfalha[j]].custoeq*PorcRep));
                    }
                    Lembrar do CIPG
                    */
                    Custoret=(Equi[posfalha[possel]].tempodereparo)*PorcRep;
                    CustoET=sal+Custoret;
                    ProdutoCEXT=CustoET*Tempo;
                    ProdutoCSXT=CustoSLA*Tempo;
                    Produtoaluguel=CustoAlugueletapa*Tempo;
                    for(int j=0;j<Num_i;j++){
                        if(Equi[j].ConsumoEnergia>0){
                            CustoEnergiadesp=CustoEnergiadesp+((Equi[j].ConsumoEnergia/1000.0)*PrKw*Tempo);
                            if(Equi[j].tipo==14||Equi[j].tipo==15){
                                CustoEnergiadespBack=CustoEnergiadespBack+((P_low/1000.0)*PrKw*Tempo);
                                //Custoextaluguel=0;
                            }
                            if(Equi[j].tipo==12){
                                int Numant;
                                if(numEq[12]>0){
                                    Numant=numEq[11];
                                }else{
                                    Numant=numEq[12];
                                }
                                CustoEnergiadespBack=CustoEnergiadespBack+(((P_high/1000.0)+Math.floor(Numant/16)*PSFP)*PrKw*Tempo);
                            }
                        }
                    }
                }
                
                //SomProd
                SomProdS=SomProdS+ProdutoCSXT;
                SomProdE=SomProdE+ProdutoCEXT;
                SomaAluguel=SomaAluguel+Produtoaluguel;
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
                /*
                if(numit==numittt){
                    FimWhile=1;
                }
                */
                //System.out.println(Arrays.toString(Eq2));
                if(Eq2[0]>=TEFOnu){
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
        double MediaCustoSLA;
        double SomaSLA=0;
        //System.out.println("custo aluguel"+custoaluguelexe);
        
        custoaluguelexe=custoaluguelexe/Num_Tentativas;
        System.out.println("custo aluguel execd"+custoaluguelexe);
        for(int i=0;i<Num_Tentativas;i++){
            SomaSLA=SomaSLA+custoits[i];
        }
        MediaCustoSLA=SomaSLA/Num_Tentativas;
        double MediaCustoEQT;
        double SomaEQT=0;
        for(int i=0;i<Num_Tentativas;i++){
            SomaEQT=SomaEQT+custoite[i];
        }
        MediaCustoEQT=SomaEQT/Num_Tentativas; 
        double VarianciaCustoSLA;
        double SomaVarianciaCustoSLA=0;
        for(int i=0;i<Num_Tentativas;i++){
            SomaVarianciaCustoSLA=SomaVarianciaCustoSLA+Math.pow((custoits[i]-MediaCustoSLA),2);
        }
        VarianciaCustoSLA=SomaVarianciaCustoSLA/Num_Tentativas;
        double VarianciaCustoEQT;
        double SomaVarianciaCustoEQT=0;
        for(int i=0;i<Num_Tentativas;i++){
            SomaVarianciaCustoEQT=SomaVarianciaCustoEQT+Math.pow((custoite[i]-MediaCustoEQT),2);
        }
        VarianciaCustoEQT=SomaVarianciaCustoEQT/Num_Tentativas;
        double DesvioPSLA=Math.sqrt(VarianciaCustoSLA);
        double DesvioPEQT=Math.sqrt(VarianciaCustoEQT);
        
        CustoEnergiadesp=CustoEnergiadesp/Num_Tentativas;
        CustoEnergiadespBack=CustoEnergiadespBack/Num_Tentativas;
        //System.out.println(Arrays.toString(custoits));
        //.println(Arrays.toString(Eq1));
        
        // CAPEX
        
        //Custo de Equipamentos
        double[] VetorCustoCapexEqui=new double[12];
        VetorCustoCapexEqui[0]=numEq[0]*preco[0];
        VetorCustoCapexEqui[1]=numEq[2]*preco[2];
        VetorCustoCapexEqui[2]=numEq[3]*preco[4];
        VetorCustoCapexEqui[3]=numEq[5]*preco[5];
        VetorCustoCapexEqui[4]=numEq[6]*preco[6];
        VetorCustoCapexEqui[5]=numEq[8]*preco[6];
        VetorCustoCapexEqui[6]=numEq[10]*preco[7];
        if(TypeMW>1){
            VetorCustoCapexEqui[7]=numEq[11]*preco[8];
            VetorCustoCapexEqui[8]=numEq[12]*preco[9];
            VetorCustoCapexEqui[9]=numEq[13]*preco[10];
            if(TypeMW>2){
                VetorCustoCapexEqui[10]=numEq[15]*preco[11];
            }
        }
        VetorCustoCapexEqui[11]=numEq[14]*preco[3];
        double CustoCompraEqui=0.0;
        for(int i=0;i<12;i++){
            CustoCompraEqui=CustoCompraEqui+VetorCustoCapexEqui[i];
        }
        double Custototinst=0;
        for(int i=0;i<numtoteq;i++){
            if(Equi[i].protecao[0]==0||Equi[i].protecao[0]>6){
                Custototinst=Custototinst+(Equi[i].instalacao+((Equi[i].distancia)/vel))*sal;
                if((Equi[i].tipo>11)&&TypeMW==1){
                Custototinst=Custototinst-(Equi[i].instalacao+((Equi[i].distancia)/vel))*sal;
                }
                if(Equi[i].tipo==12&&TypeMW==2){
                Custototinst=Custototinst-(Equi[i].instalacao+((Equi[i].distancia)/vel))*sal;
                }
            }
        }
        //Custos de Fibra
        double sumaux1=0;
        double sumaux2=0;
        for(int i=0;i<(n/2.0);i++){
            for(int j=0;j<(n/2.0);j++){
                sumaux1=sumaux1+(((n-2*i-1)/2)+((n-2*j-1)/2));
            }
        }
        for(int i=0;i<(N/2.0);i++){
            for(int j=0;j<(N/2.0);j++){
                sumaux2=sumaux2+(((N-2*i-1)/2)+((N-2*j-1)/2));
            }
        }
        double Custotfiber=(4*sumaux1*l*N*N+4*sumaux2*L)*PkmVF;
        double Custottren=(((n*n-1)*l*N*N)+((N*N-1)*L))*PkmTF;
        double CustoTfiberPFL=0;
        double CustoTfiberPDL=0;
        if(selENA==2){
            CustoTfiberPFL=(N-1)*N*L*PkmTF;
            Custotfiber=Custotfiber+((Math.pow((N/2.0),2)+0.5*N)+(Math.pow((N/2.0),2)-0.5*N))*N*L*(Math.ceil(((n*n)/(SR))))*PkmVF;
        }
        if(selEND==2){
            CustoTfiberPDL=(n-1)*n*N*N*l*PkmTF;
            Custotfiber=Custotfiber+((Math.pow((n/2.0),2)+0.5*n)+(Math.pow((n/2.0),2)-0.5*n))*n*l*N*N*PkmVF;
            
        }
        //Custos de Aluguel de Esptrum
        double Custo_aluguelEsp = 0;
        if(TypeMW==2||TypeMW==3){
            Custo_aluguelEsp=(numEq[11]+numEq[12])*2*PrSpleaYe*Ny;
        }
        double Custocapex=Custotfiber+Custottren+Custototinst+CustoCompraEqui+CustoTfiberPFL+CustoTfiberPDL+Custo_aluguelEsp;
        //custos de Energia
        for(int i=0;i<numtoteq;i++){
            Custo_Energia=Custo_Energia+Equi[i].ConsumoEnergia;
        }
        //System.out.println(Custo_Energia/1000.0);
        Custo_Energia=(Custo_Energia/1000.0)*Ny*365*24*PrKw;
        if(Custo_Energia>0){
            Custo_Energia=Custo_Energia-CustoEnergiadesp;
        }
        double CustoBack= ((numEq[11]+numEq[12])*(P_low/1000.0)+numEq[15]*(P_high/1000.0)+Math.floor((numEq[11]+numEq[12])/16)*PSFP)*Ny*365*24*PrKw;;
        CustoBack=CustoBack-CustoEnergiadespBack;
        if(TypeMW==3){
            Custo_Energia=Custo_Energia+CustoBack;
        }
        //Custos de Aluguel Gb/s
        double Custo_aluguel=0;
        if(TypeMW==1){
            Custo_aluguel=((Nummedclifal*100)/1024.0)*PrGbsh*24*365*Ny;
            System.out.println("aluguel paad"+Custo_aluguel);
            Custo_aluguel=custoaluguelexe+Custo_aluguel;
        }
        ///Resultado
        String tiprot1="";
        String tiprot2="";
        String tiprot3="";
        String tiprot4="";
        String tiprot5="";
        String timw="";
        if(selCO==1){
            tiprot1="OLT chassi+portas";
        }else if(selCO==2){
            tiprot1="OLT protegido";
        }
        if(selENA==1){
            tiprot2="alimentação simples";
        }else if(selCO==2){
            tiprot2="alimentação dupla";
        }
        if(selRN==1){
            tiprot3="splitter padrão";
        }else if(selRN==2){
            tiprot3="AWG";
        }
        if(selEND==1){
            tiprot4="distribuição simples";
        }else if(selEND==2){
            tiprot4="distribuição dupla";
        }
        if(selRES==1){
            tiprot5="ONU";
        }else if(selRES==2){
            tiprot5="ONU com MW";
        }
        if(TypeMW==0){
            timw="Sem MW";
        }else if(TypeMW==1){
            timw="MW leasing";
        }else if(TypeMW==2){
            timw="MW upgrade";
        }else if(TypeMW==3){
            timw="MW propietario";
        }
        ///Resultado
        System.out.println( tiprot1 +" , "+ tiprot2 +" , "+tiprot3 +" , "+tiprot4 +" , " +tiprot5 +" , " + timw + " : ");
        System.out.println( "Custo de Penalidade:  " + MediaCustoSLA);
        System.out.println( "Custo de Reparo:  " + MediaCustoEQT);
        System.out.println( "Custo de Inst de Eq:  " + Custototinst);
        System.out.println( "Custo de Compra de Eq:  " + CustoCompraEqui);
        System.out.println( "Custos Relacionados a Fibra:  " + (Custocapex-CustoCompraEqui-Custototinst));
        System.out.println( "Custo de Aluguel de Espectro:  " + Custo_aluguelEsp);
        System.out.println( "Capex:  " + Custocapex);
        System.out.println( "Custo de Aluguel:  " + Custo_aluguel);
        System.out.println( "Custo de Energia:  " + Custo_Energia);
        System.out.println( "Média de falhas de equipamentos:  "+"{Onu,FDS,Splitters,Olt ports,Rn_Chassi,FFS,Olt_Chassi,Olt_Switch,PFFS,Olt_Switch,PDFS,GES,Macro,Antena,Micro,Pico}");
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
    private double ObterAluguel(int numcli,int Num_i,double[][] Matriz_Parametro,double PorcentagemCom,int[][] CenarioMWCli,int[][]Cenaricli,int nummed){
        int NF=0;
        int[] VF=new int[numcli];
        int[][] VFfin=new int[Num_i][numcli];
        int[] VFred=new int[2];
        int[] VFaux=new int[numcli];
        int[] VPAred=new int[2];
        int[] VPAaux=new int[numcli];
        int[] VPA=new int[numcli];
        int[] VPDred=new int[2];
        int[] VPDaux=new int[numcli];
        int[] VPD=new int[numcli];
        int[] VPOLT=new int[numcli];
        int[] VPOLTaux=new int[numcli];
        int[] VPOLTaux2=new int[numcli];
        //especificação da proteção
        for(int i=0;i<Num_i;i++){
            if(TypeMW==1){
                if((Matriz_Parametro[4][i]==2)||(Matriz_Parametro[4][i]==3)){
                    VPAred[0]=(int) Matriz_Parametro[5][i];
                    VPAred[1]=(int) Matriz_Parametro[6][i];       
                    VPAaux=ObterVetorCompleto(VPAred);
                    if(Matriz_Parametro[4][i]==2){
                        VPAaux=ObterVetorPrt2(numcli,VPAaux);
                    }else if(Matriz_Parametro[4][i]==3){
                        VPAaux=ObterVetorPrt(numcli,VPAaux);
                    }
                    for(int j=0;j<numcli;j++){
                        if((VPA[j]+VPAaux[j])>0){
                            VPA[j]=1;
                        }
                    }
                    VPAaux=new int[numcli];
                }else if((Matriz_Parametro[4][i]==5)||(Matriz_Parametro[4][i]==6)){
                    VPDred[0]=(int) Matriz_Parametro[5][i];
                    VPDred[1]=(int) Matriz_Parametro[6][i];       
                    VPDaux=ObterVetorCompleto(VPDred);
                    if(Matriz_Parametro[4][i]==5){
                        VPDaux=ObterVetorPrt2(numcli,VPDaux);
                    }else if(Matriz_Parametro[4][i]==6){
                        VPDaux=ObterVetorPrt(numcli,VPDaux);
                    }
                    for(int j=0;j<numcli;j++){
                        if((VPD[j]+VPDaux[j])>0){
                            VPD[j]=1;
                        }
                    }
                    VPDaux=new int[numcli];
                }else if(Matriz_Parametro[4][i]==8){
                    int [] MP={0,0};
                    MP[0]=(int) Matriz_Parametro[5][i];
                    MP[1]=(int) Matriz_Parametro[6][i];
                    VPOLTaux=ObterVetorCompleto(MP);
                    for(int j=0;j<numcli;j++){
                        if((VPOLT[j]+VPOLTaux[j])>0){
                            VPOLT[j]=1;
                        }
                    }

                }
            }
        }
	//daqui começa
	for(int i=0;i<Num_i;i++){
            if(TypeMW==1){
		VFred[0]=(int) Matriz_Parametro[0][i];
		VFred[1]=(int) Matriz_Parametro[1][i];       
		VFaux=ObterVetorCompleto(VFred);
		int a=(int) Matriz_Parametro[8][i];
		if(a==0||a==1||a==2||a==3||a==4||a==5||a==6||a==7||a==9||a==11||a==16||a==17){
			for(int j=0;j<numcli;j++){
                	    VFfin[i][j]=VFaux[j];
               		}
		}
		if((a==3||a==6||a==17)&&selCO==2){
                        if(VFred[0]==0&&VFred[0]==0){
                            for(int j=0;j<numcli;j++){
                                VFfin[i][j]=0;
                            }
			}else{
                            for(int j=0;j<numcli;j++){
                                    if((VFfin[i][j]+VPOLT[j])==2){
                                            VFfin[i][j]=1;
                                    }else{
                                    VFfin[i][j]=0;
                                }
                            }	
                    }
                }
		if((a==5)&&selENA==2){
			for(int j=0;j<numcli;j++){
                		if((VFfin[i][j]+VPA[j])==2){
					VFfin[i][j]=1;
				}else{
                                    VFfin[i][j]=0;
                                }
               		}	
		}
		if((a==1)&&selEND==2){
                    System.out.println(Arrays.toString(VPD));
			for(int j=0;j<numcli;j++){
                		if((VFfin[i][j]+VPD[j])==2){
					VFfin[i][j]=1;
				}else{
                                    VFfin[i][j]=0;
                                }
               		}	
		}
            }
	}
	for(int i=0;i<Num_i;i++){
		for(int j=0;j<numcli;j++){
            		VF[j]=VF[j]+VFfin[i][j];
			if(VF[j]>1){
				VF[j]=1;	
			}
            	}
        }
	for(int i=0;i<numcli;i++){
            if(VF[i]==1){
                NF++;
            }
        }
        int NumExe=NF-nummed;
        double Custo_aluguel=0;
        if(NF>nummed){
            Custo_aluguel=Math.pow((NF-nummed),paraluext)*PrGbsh;
            //Custo_aluguel=(((NF-nummed)*100)/1024)*PrGbsh*((228*Math.exp(0.001*(NF-nummed)))-220);
            //System.out.println("custo exp:  "+Custo_aluguel);
            //System.out.println("custo lin:  "+((((NF-nummed)*100)/1024)*PrGbsh*(NF-nummed)));
        }
        //System.out.println("custo tot:  "+Custo_aluguel);
        return Custo_aluguel;
    }
    private double ObterCusto(int numcli,int Num_i,double[][] Matriz_Parametro,double PorcentagemCom,int[][] CenarioMWCli,int[][]Cenaricli){
        int NF=0;
        int[] VF=new int[numcli];
        int[][] VFfin=new int[Num_i][numcli];
        int[] VFred=new int[2];
        int[] VFaux=new int[numcli];
        int[] VPAred=new int[2];
        int[] VPAaux=new int[numcli];
        int[] VPA=new int[numcli];
        int[] VPDred=new int[2];
        int[] VPDaux=new int[numcli];
        int[] VPD=new int[numcli];
        int[] VPMW=new int[numcli];
        int[] VPMWaux=new int[numcli];
        int[] VPOLT=new int[numcli];
        int[] VPOLTaux=new int[numcli];
        int[] VPOLTaux2=new int[numcli];
        //especificação da proteção
        for(int i=0;i<Num_i;i++){
            if(Matriz_Parametro[7][i]==1.0){
                if((Matriz_Parametro[4][i]==2)||(Matriz_Parametro[4][i]==3)){
                    VPAred[0]=(int) Matriz_Parametro[5][i];
                    VPAred[1]=(int) Matriz_Parametro[6][i];     
                    VPAaux=ObterVetorCompleto(VPAred);
                    if(Matriz_Parametro[4][i]==2){
                        VPAaux=ObterVetorPrt2(numcli,VPAaux);
                    }else if(Matriz_Parametro[4][i]==3){
                        VPAaux=ObterVetorPrt(numcli,VPAaux);
                    }
                    for(int j=0;j<numcli;j++){
                        if((VPA[j]+VPAaux[j])>0){
                            VPA[j]=1;
                        }
                    }
                    VPAaux=new int[numcli];
                }else if((Matriz_Parametro[4][i]==5)||(Matriz_Parametro[4][i]==6)){
                    VPDred[0]=(int) Matriz_Parametro[5][i];
                    VPDred[1]=(int) Matriz_Parametro[6][i];       
                    VPDaux=ObterVetorCompleto(VPDred);
                    if(Matriz_Parametro[4][i]==5){
                        VPDaux=ObterVetorPrt2(numcli,VPDaux);
                    }else if(Matriz_Parametro[4][i]==6){
                        VPDaux=ObterVetorPrt(numcli,VPDaux);
                    }
                    for(int j=0;j<numcli;j++){
                        if((VPD[j]+VPDaux[j])>0){
                            VPD[j]=1;
                        }
                    }
                    VPDaux=new int[numcli];
                }else if((Matriz_Parametro[4][i]==7)&&(TypeMW>1)){
                    if(Matriz_Parametro[8][i]==12){
                        for(int j=0;j<numcli;j++){
                            VPMWaux[j]=1;
                        }
                    }else{
                        int MP1=(int) Matriz_Parametro[5][i];
                        int MP2=(int) Matriz_Parametro[6][i];
                        for(int t=MP1;t<=MP2;t++){
                            for(int j=0;j<n*N;j++){
                                for(int k=0;k<n*N;k++){
                                    if((CenarioMWCli[j][k]>=Matriz_Parametro[5][i])&&(CenarioMWCli[j][k]<=Matriz_Parametro[6][i])){
                                        VPMWaux[Cenaricli[j][k]]=1;
                                    }
                                }
                            }
                        }
                    }
                    for(int j=0;j<numcli;j++){
                        if((VPMW[j]+VPMWaux[j])>0){
                            VPMW[j]=1;
                        }
                    }
                }else if(Matriz_Parametro[4][i]==8){
                    int [] MP={0,0};
                    MP[0]=(int) Matriz_Parametro[5][i];
                    MP[1]=(int) Matriz_Parametro[6][i];
                    VPOLTaux=ObterVetorCompleto(MP);
                    for(int j=0;j<numcli;j++){
                        if((VPOLT[j]+VPOLTaux[j])>0){
                            VPOLT[j]=1;
                        }
                    }

                }
            }
        }
	//daqui começa
	for(int i=0;i<Num_i;i++){
            if(Matriz_Parametro[7][i]==1.0){
		VFred[0]=(int) Matriz_Parametro[0][i];
		VFred[1]=(int) Matriz_Parametro[1][i];       
		VFaux=ObterVetorCompleto(VFred);
		int a=(int) Matriz_Parametro[8][i];
		if(a==0||a==1||a==2||a==3||a==4||a==5||a==6||a==7||a==9||a==11||a==16||a==17){
			for(int j=0;j<numcli;j++){
                	    VFfin[i][j]=VFaux[j];
               		}
		}
		if((a==3||a==6||a==17)&&selCO==2){
                        if(VFred[0]==0&&VFred[0]==0){
                            for(int j=0;j<numcli;j++){
                                VFfin[i][j]=0;
                            }
			}else{
                            for(int j=0;j<numcli;j++){
                                    if((VFfin[i][j]+VPOLT[j])==2){
                                            VFfin[i][j]=1;
                                    }else{
                                    VFfin[i][j]=0;
                                }
                            }	
                    }
                }
		if((a==5)&&selENA==2){
			for(int j=0;j<numcli;j++){
                		if((VFfin[i][j]+VPA[j])==2){
					VFfin[i][j]=1;
				}else{
                                    VFfin[i][j]=0;
                                }
               		}	
		}
		if((a==1)&&selEND==2){
                    //System.out.println(Arrays.toString(VPD));
			for(int j=0;j<numcli;j++){
                		if((VFfin[i][j]+VPD[j])==2){
					VFfin[i][j]=1;
				}else{
                                    VFfin[i][j]=0;
                                }
               		}	
		}
		if((a!=0)&&(TypeMW>0)){
			for(int j=0;j<numcli;j++){
                		if((VFfin[i][j]+VPMW[j])==2){
                                    VFfin[i][j]=1;
				}else{
                                    VFfin[i][j]=0;
                                }
               		}	
		}
            }
	}
        double NFC=0;
        double NFR=0;
        double CustoSLA=0;
	for(int i=0;i<Num_i;i++){
		for(int j=0;j<numcli;j++){
            		VF[j]=VF[j]+VFfin[i][j];
			if(VF[j]>1){
				VF[j]=1;	
			}
            	}
        }
	for(int i=0;i<numcli;i++){
            if(VF[i]==1){
                NF++;
            }
        }
        //System.out.println("NF: "+NF);
        NFC=NF*PorcentagemCom;
        NFR=NF-NFC;
        CustoSLA=Math.pow(NFC,alpha)*SLAC+Math.pow(NFR,alpha)*SLAR;
        //tem.out.println("custo tot:  "+CustoSLA);
        return CustoSLA;
    }
}


