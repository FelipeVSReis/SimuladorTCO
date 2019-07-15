package SimuladorTCO;


import java.util.Arrays;
import java.util.Random;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author felip
 */
public class Hexagono {
    double l = 0.03;
    int n = 10;
    int N = 10;
    int SR = 16; //Limitado a 16,32 ou 64 
    
    public Hexagono(double raio,int parmetodo){
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
        /*
        //Para debug
        int [] vetorDeImpressao = new int[cobertura[0].length];
        for (int i=0; i<cobertura.length;i++){
            for (int j=0; j<cobertura[0].length; j++){
            vetorDeImpressao[j] = cobertura [i][j];
        }
        System.out.println(Arrays.toString(vetorDeImpressao));
    }
        */
        int teste = contaCelulas(cobertura);
        
        int [][]matrizDeCobertura = matrizDePosicoes(cobertura, teste, 16);
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
        
        matSimplificada = MatSimplificada(matSimplificada, SR);
        
        int cont2=0;
        for(int i=0;i<matSimplificada.length;i++){
            for (int j=0; j<matSimplificada[0].length;j++){
                if(matSimplificada[i][j]>0){
                    matDePos[cont2][2] = matSimplificada[i][j];
                    cont2++;
                }
            }
        }
        
        int [] vetorDeImpressao2 = new int[matDePos[0].length];
        for (int i=0; i<matDePos.length;i++){
            for (int j=0; j<matDePos[0].length; j++){
            vetorDeImpressao2[j] = matDePos [i][j];
            }
        System.out.println(Arrays.toString(vetorDeImpressao2));
        }
        System.out.println("termina aqui2");
        
        return matDePos;
    }
     
    private int[][] MatSimplificada (int[][] matSimplificada,int SR){
        int delH, delV, vIndex = 1, vIndexReset=1, auxH=0;
        if (SR ==16){
        delH=4;
        delV=4;
        } else if(SR==32){
        delH=4;
        delV=8;
        } else if (SR == 64){
        delH=8;
        delV=8;
        } else if (SR == 8){
        delH = 2;
        delV = 4;
        } else {
        delH = 0;
        delV = 0;
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

