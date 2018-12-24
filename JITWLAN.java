// writer : Kimin Lee
// JIT-WLAN algorithm in the paper 'http://ieeexplore.ieee.org/document/7524451/'

import java.sql.Date;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.ArrayList;
import java.util.concurrent.ConcurrentSkipListSet;
import java.util.Set;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;

public class JITWLAN {
	
	// variables 
	private int Num_user, Num_AP, Num_Ch;
	private double trade_off_H;
	private double [][] Weight; // W_{ij} = (queue_j) X (archievable rate of j from i)
	private double [] Cost_AP; // cost to turn on AP j
	private int [][] Overlap_G; // overlap graph
	private int [] operating_channel; // operating channel of AP j

	/* generator
	 * input a = # of users, b = # of APs, c = # of non overlapping channels and d = tradeoff parameter H
	 * */
	public JITWLAN(int N_uer, int N_AP, int N_CH, double t_H){
		Num_user = N_uer;
		Num_AP = N_AP;
		Num_Ch = N_CH;
		trade_off_H = t_H;
		Cost_AP = new double[Num_AP];
		Weight = new double[Num_AP][Num_user];
		Overlap_G = new int[Num_AP][Num_AP];
		operating_channel = new int[Num_AP];
	}
	
	/* 
	 * reset the variables
	 * */
	public void initialize(int N_uer, int N_AP, int N_CH, double t_H){
		Num_user = N_uer;
		Num_AP = N_AP;
		Num_Ch = N_CH;
		trade_off_H = t_H;
		Cost_AP = new double[Num_AP];
		Weight = new double[Num_AP][Num_user];
		Overlap_G = new int[Num_AP][Num_AP];
		operating_channel = new int[Num_AP];
	}
	
    public boolean check_weight(){
        boolean result = false;
        for (int i =0;i<Num_AP;i++){
               for (int j =0; j<Num_user;j++){
                      if(Weight[i][j]>0){
                            result=true;
                            break;
                      }
               }
        }
        return result;
    }
    
	// get operating channel of APs
	public int [] get_oper_channel(){
		return operating_channel;
	}
	
	// get the value of Num_user
	public int getNum_user(){
		return Num_user;
	}
	
	// get the value of Num_AP
	public int getNum_AP(){
		return Num_AP;
	}
	
	// get the value of Num_Ch
	public int getNum_Ch(){
		return Num_Ch;	
	}
	
	// get the value of Num_Ch
	public double get_trade_off_H(){
		return trade_off_H;	
	}
	
	// get the value of weight[a][b]
	public double getWeight(int a, int b){
		return Weight[a][b];
	}
	
	// set the value of weight[a][b]= c
	public void setWeight(int a, int b, double c){
		Weight[a][b] = c;
	}
	
	// get the value of Overlap_G[a][b]
	public int getOverlap_G(int a, int b){
		return Overlap_G[a][b];
	}
	
	// set the value of Overlap_G[a][b]= c
	public void setOverlap_G(int a, int b, int c){
		Overlap_G[a][b] = c;
	}
	
	// get the value of Cost_AP[a]
	public double getCostAP(int a){
		return Cost_AP[a];
	}
	
	// set the value of Cost_AP[a]= b
	public void setCost(int a, double b){
		Cost_AP[a] = trade_off_H*b;
	}
	
	// return the real index of virtual AP 
	public int get_real_index(int index, int num_Ch){
		return (int) Math.ceil(index/num_Ch);
	}
	
	// return a operating channel implied by virtual AP
	public int get_channel_index (int index, int num_Ch){
		int return_value =(index+1)%num_Ch;
		if(return_value==0){
			return_value = num_Ch;
		}	
		return return_value;
	}
	
	// generate a virtual interference graph
	public int[][] virtual_interference_graph(int [][]O_G, int num_Ch, int num_AP){
		int V_Num_AP = num_AP*num_Ch;
		int [][] return_AP = new int[V_Num_AP][V_Num_AP];
		
		for(int i =0;i<V_Num_AP; i++){
			int real_index = get_real_index(i,num_Ch);
			int channel_index = get_channel_index(i,num_Ch);
			
			for (int j = 0; j<V_Num_AP;j++){
				// interference with same APs
				if((real_index*num_Ch<=j)&&(j<real_index*num_Ch+num_Ch)){
					return_AP[i][j] =1;
				}
				// interference with other APs
				if(O_G[real_index][get_real_index(j,num_Ch)]==1){
					if(channel_index==get_channel_index(j,num_Ch)){
						return_AP[i][j] =1;
					}
				}
			}
			return_AP[i][i] =0;
		}
		return return_AP;
	}
	
	// calculate the d_max (maximum # of users which are served by one AP)
	public int Cal_d_max (int num_AP, int num_user, double [][] weight){
		int d_max = 0;
		for (int i =0;i<num_AP;i++){
			int temp_max = 0;
			for (int j =0; j<num_user;j++){
				if(weight[i][j]>0){
					temp_max +=1;
				}
			}
			if(temp_max>d_max){
				d_max =temp_max;
			}
		}
		return d_max;
	}
	
	// find an optimal solution of SLP using CPLEX
	// return the solution from CPLEX
	public double [] solveLP(){
		
		double [] LP_sol = new double [Num_user*Num_AP*Num_Ch+Num_AP*Num_Ch]; // solution of LP
		try {
			int d_max = Cal_d_max(Num_AP,Num_user,Weight);
			int[][] V_I_AP = virtual_interference_graph(Overlap_G,Num_Ch,Num_AP);
			IloCplex cplex = new IloCplex();
			
			// association variable & on_off variable
			// z11 z21 .. zU1 z12 z22 .. zU2 .. z1N z2N .. zUN y11 y12 ..y1c y21 y22 .. y2c .. yN1 yN2 .. yNc
			IloNumVar [] z = cplex.numVarArray(Num_user*Num_AP*Num_Ch + Num_AP*Num_Ch, 0.0, 1);
		
			// set the coefficient for z
			double [] W = new double[Num_user*Num_AP*Num_Ch+Num_AP*Num_Ch];		
			
			for (int j =0;j<Num_AP;j++){
				for (int c = 0; c<Num_Ch;c++){
					W[Num_user*Num_AP*Num_Ch + j*Num_Ch+c] = -Cost_AP[j];
					for(int i = 0;i<Num_user;i++){
						W[j*Num_Ch*Num_user+Num_user*c+i]= Weight[j][i];
					}
				}
			}
			
			// objective fn
			IloLinearNumExpr objective = cplex.scalProd(z, W);
			
			// max obj
			cplex.addMaximize(objective);
			
			// constraints 
			// single assoication & cover all users
			if(d_max ==0){
				return LP_sol;
			}
			
			double epsilon = (double) 1/(d_max);
			for (int i =0;i<Num_user; i++){
				double [] coeff = new double[Num_user*Num_AP*Num_Ch+ Num_AP*Num_Ch];
				for (int j =0;j<Num_AP*Num_Ch;j++){
					coeff[i+j*Num_user] =1;
				}
				// sum_{j} zij <1(lambda)
				cplex.addLe(cplex.scalProd(z, coeff),1);
				// sum_{j} zij > epsilon
				System.out.println("epsilon : " + epsilon + "  d_max : " + d_max);
				cplex.addGe(cplex.scalProd(z, coeff),epsilon);
				 
			}
			
			//configuration 
			for (int j=0;j<Num_AP*Num_Ch;j++){
				double [] coeff2 = new double[Num_user*Num_AP*Num_Ch+ Num_AP*Num_Ch];
				for (int i=0; i<Num_user;i++){
					coeff2[j*Num_user +i]=1;
				}
				coeff2[Num_user*Num_AP*Num_Ch+j] =-1;
				//sum_{i} zij < sum_{c} yjc
				cplex.addLe(cplex.scalProd(z, coeff2),0);
			}
			
			//interference
			for (int j=0;j<Num_AP*Num_Ch;j++){
				for (int i =0; i<Num_AP*Num_Ch;i++){
					if(V_I_AP[j][i]==1){
						double [] coeff3 = new double[Num_user*Num_AP*Num_Ch+ Num_AP*Num_Ch];
						coeff3[Num_user*Num_AP*Num_Ch+i] =1;
						coeff3[Num_user*Num_AP*Num_Ch+j] =1;
						cplex.addLe(cplex.scalProd(z, coeff3),1);
					}						
				}
			}
			 
			// solve
			if(cplex.solve()){
				System.out.println("model solved");
				System.out.println("obj = "+cplex.getObjValue());
				LP_sol = cplex.getValues(z);
			}else{
				System.out.println("model not solved");
			}
			
		} catch (IloException exc) {
			exc.printStackTrace();
		}
		return LP_sol;
	}
	
	// find the cost value of given input
	public double Cost_function(int [][] solution){
		double return_value =0;
		int [] count_AP_user = new int [Num_AP]; // the # of user per each APs
		// count the # of users per each APs
		for (int i =0; i<Num_AP; i++){
			for(int j=0; j<Num_user; j++){
				if(solution[i][j]==1){
					count_AP_user[i] += 1;
				}					
			}
			// add cost
			if(count_AP_user[i] !=0){
				return_value -= Cost_AP[i]; // - cost to turn on AP
				for(int j=0; j<Num_user;j++){
					if(solution[i][j]==1){
						return_value += (double) Weight[i][j]/count_AP_user[i];
					}					
				}
			}
		}
		return return_value;
	}
	
	// calculate the degree of AP 
	public int [] Degree_LP_graph(int[][] V_I_AP, int [] num_user_per_AP){
		
		int V_num_AP = Num_AP*Num_Ch;
		int [] degree = new int [V_num_AP];
		
		for (int i =0; i<V_num_AP; i++){
			if(num_user_per_AP[i]>0){ // if AP i is selected
				for(int j =0; j<V_num_AP; j++){
					if((V_I_AP[i][j]==1)&&(num_user_per_AP[j]>0)){
						degree[i] +=1;
					}
				}
			} 
			// otherwise, degree of AP i = 0
		}
		return degree;
	}
	
	// calculate the theta bar
	public double [] Generate_theta_bar(double [] theta_j, int [] degree){
		int V_num_AP = Num_AP*Num_Ch;
		double [] theta_bar = new double [V_num_AP];
		
		for (int i =0; i<V_num_AP; i++){
			if(theta_j[i]>0){
				theta_bar[i] = (double) theta_j[i]/(degree[i]+1);
			}
		}
		
		return theta_bar;
	}
	
	public int sum_array(int [] target){
		int value=0;
		for(int i=0;i<target.length;i++){
			value += target[i];
		}
		return value;
	}
	
	// select the on AP using degree based greedy algorithm and user association
	public int [][] Round_max(double[] LP_sol){
		int V_num_AP = Num_AP*Num_Ch; // total # of virtual AP
		double [][] Weight_matrix = Weight.clone(); // modified wieght_matrix
		int [] num_user_per_AP = new int[V_num_AP]; // # of users per AP
		int [] V_on = new int[V_num_AP]; // set of on APs
		int [] V_off = new int[V_num_AP]; // set of off APs
		int[][] V_I_AP = virtual_interference_graph(Overlap_G,Num_Ch,Num_AP); // virtual interference
		
		double [] theta_j = new double [V_num_AP]; // theta from optimal
		int [] V_star = new int[V_num_AP];
		 
		// calculate theta
		for (int i =0; i<V_num_AP;i++){
			int real_index = get_real_index(i,Num_Ch);
			int [][]temp_sol  = new int [Num_AP][Num_user];
			for (int j=0;j <Num_user;j++){
				if(LP_sol[i*Num_user+j]>0){
					temp_sol[real_index][j] =1;
					num_user_per_AP[i]+=1;
					V_star[i]=1;
				}				
			}
			theta_j[i] = Cost_function(temp_sol);
		}
		
		while(true){
			int sum = sum_array(V_star);
            if(sum==0){
                   System.out.println("AP selection is completed");
                   break;
            }
			
			int [] degree = Degree_LP_graph(V_I_AP,num_user_per_AP);
			// calculate theta bar
			double [] theta_bar_j = Generate_theta_bar(theta_j,degree);
			
			// find optimal node
			double max_val =Double.NEGATIVE_INFINITY;
			int max_index =0;
			
			for (int i =0; i<V_num_AP;i++){
				if(V_star[i]==1){
					if(Double.compare(theta_bar_j[i], max_val)==1){ // maximum case
						max_val = theta_bar_j[i];
						max_index = i;
					}
				}
			}
			
			// turn on selected AP
			V_on[max_index] =1;
			V_star[max_index]=0;
			
			//turn off neighbor
			for (int i=0; i<V_num_AP;i++){
				if(V_I_AP[max_index][i]==1){ //neighbor node
					// turn off 
					V_off[i] =1;
					V_star[i]=0;
					// set z =0
					for(int j=0;j<Num_user;j++){
						LP_sol[i*Num_user+j] =0;
					}
				}
			}
		}
		
		int [] U_on =new int [Num_user]; // U_on: the set of users who have positive z
		int [][] fix_sol = new int[Num_AP][Num_user]; // final solution
		double [][] V_cost_AP = new double [Num_user][V_num_AP]; // virtual cost when user i select virtual AP j
		
		for (int i =0; i<Num_user;i++){
			for (int j =0;j<V_num_AP;j++){
				
				if(LP_sol[j*Num_user+i]>0){ // if it has positive Z*
					U_on[i]=1;
					V_cost_AP[i][j] = (double) Weight[get_real_index(j,Num_Ch)][i]/num_user_per_AP[j]; //sigma
				}
			}
			
			// assign U_on
			if(U_on[i]==1){
				// find optimal node
				double max_val =0;
				int max_index =0;
				
				for (int j =0;j<V_num_AP;j++){
					if(Double.compare(V_cost_AP[i][j], max_val)==1){ // maximum case
						max_val = V_cost_AP[i][j];
						max_index = j;
					}
				}
				//assign user
				fix_sol[get_real_index(max_index,Num_Ch)][i] =1;
				operating_channel[get_real_index(max_index,Num_Ch)] =get_channel_index (max_index,Num_Ch);
			}
		}
		
		//assign U_off
		for (int i =0; i<Num_user;i++){
			// assign U_off
			if(U_on[i]==0){
				// calculate the delta 
				for (int j =0;j<V_num_AP;j++){
					int real_index_j =get_real_index(j,Num_Ch);
					
					if((V_off[j]==0)&&(Weight[real_index_j][i]>0)){
						int [][] temp_sol = new int [Num_AP][Num_user];
						// copy array
						for(int k=0;k<Num_AP;k++){
							for(int kk=0;kk<Num_user;kk++){
								temp_sol[k][kk] = fix_sol[k][kk];
							}
						}
						temp_sol[real_index_j][i] = 1;
						V_cost_AP[i][j] = Cost_function(temp_sol); // delta
					}
				}
				
				// find max node
				double max_val =Double.NEGATIVE_INFINITY;
				int max_index =0;
				
				for (int j =0;j<V_num_AP;j++){
					if(V_off[j]==0){
						if(Double.compare(V_cost_AP[i][j], max_val)==1){ // maximum case
							max_val = V_cost_AP[i][j];
							max_index = j;
						}
					}
				}
				
				int real_index_j =get_real_index(max_index,Num_Ch);
				
				//assign user
				fix_sol[real_index_j][i] =1;
				operating_channel[real_index_j] =get_channel_index (max_index,Num_Ch);
				
				//if it is not selected AP
				if(V_on[max_index]==0){
					// turn on
					V_on[max_index] =1;
					//turn off neighbor
					
					for (int j =0;j<V_num_AP;j++){
						if(V_I_AP[max_index][j]==1){ //neighbor node
							V_off[j] =1;
						}
					}
				}
			}
		}
		
		return fix_sol;
	}
	
	public int [][] Round_min(double[] LP_sol){

		int V_num_AP = Num_AP*Num_Ch; // total # of virtual AP
		double [][] Weight_matrix = new double [Num_AP][Num_user]; // modified wieght_matrix
		int [] user_selected = new int[Num_user]; // flag for user
		int [][] fix_sol = new int[Num_AP][Num_user]; // final solution
		int [] V_on = new int[V_num_AP]; // set of on APs
		int [] V_off = new int[V_num_AP]; // set of off APs
		int[][] V_I_AP = virtual_interference_graph(Overlap_G,Num_Ch,Num_AP); // virtual interference
		
		for (int i =0; i<Num_user;i++){
			for (int j =0;j<Num_AP;j++){
				Weight_matrix[j][i] = Weight[j][i];
			}
		}
		
		// remove unavailable users
		for (int i =0; i<Num_user;i++){
			double sum_weight = 0;
			for (int j =0;j<Num_AP;j++){
				sum_weight += Weight_matrix[j][i];
			}
			if(sum_weight==0){
				user_selected[i] = 1;
			}
		}
		
		//check the weight
		if(sum_array(user_selected)==Num_user){
			System.out.println("every user has 0 weight");
			return fix_sol;
		}

		// select APs
		for (int j =0; j<V_num_AP;j++){
			if(LP_sol[j]>0.5){
				int real_index = get_real_index(j,Num_Ch);
				V_on[j] = 1;
				//turn off neighbor
				for (int i=0; i<V_num_AP;i++){
					if(V_I_AP[j][i]==1){ //neighbor node
						// turn off 
						V_off[i] =1;
					}
				}
				
				//user modification
				for (int i =0; i<Num_user;i++){
					//covered user
					if(Weight_matrix[real_index][i]>0){
						user_selected[i] = 1;
						//update Weight
						for (int k =0; k<Num_AP;k++){
							Weight_matrix[k][i] =0;
						}
					}
				}
			}
		} 
		
		//check the weight
		if(sum_array(user_selected)!=Num_user){
			//greedy selection
			while(true){
				double [] metric = new double [V_num_AP];
				
				for (int j =0; j<V_num_AP;j++){
					
					// for available ap
					if((V_on[j]==0)&&(V_off[j]==0)){
						
						int real_index = get_real_index(j,Num_Ch);
						int available_user = 0;
						
						//cal available users
						for (int i =0; i<Num_user;i++){
							if(Weight_matrix[real_index][i]>0){
								available_user+=1;
							}
						}
						
						if(available_user!=0){
							metric[j] =Cost_AP[real_index]/(available_user);
						}else{
							metric[j] =Double.NEGATIVE_INFINITY;
						}
						
					} else{
						metric[j] =Double.NEGATIVE_INFINITY;
					}
				}
				
				//select max
				int max_index = 0;
				double max_val = Double.NEGATIVE_INFINITY;
				for (int i =0; i<V_num_AP;i++){						
					if(Double.compare(metric[i], max_val)==1){ // maximum case
						max_val = metric[i];
						max_index = i;
					}
				}
				
				int real_index = get_real_index(max_index,Num_Ch);
				
//					System.out.println("selected v AP" +max_index);
				V_on[max_index] = 1;
				
				//turn off neighbor
				for (int i=0; i<V_num_AP;i++){
					if(V_I_AP[max_index][i]==1){ //neighbor node
						// turn off 
						V_off[i] =1;
					}
				}
				
				//user modification
				for (int i =0; i<Num_user;i++){
					//covered user
					if(Weight_matrix[real_index][i]>0){
						user_selected[i] = 1;
						//update Weight
						for (int k =0; k<Num_AP;k++){
							Weight_matrix[k][i] =0;
						}
					}
				}
				
				if(sum_array(user_selected)==Num_user){
					System.out.println("cover_every_user_by_greedy");
					break;
				}
			}
			
		}else{
			System.out.println("cover_every_user_by_SCLP");
		}
		
		
		for(int i=0; i<V_num_AP;i++){
			System.out.println("V_on["+i+"] = "+V_on[i]);
		}

		double [] cover_lp = set_cover_rounding(V_on);
		
		for(int i=0; i<Num_user*Num_AP*Num_Ch+Num_AP*Num_Ch;i++){
			System.out.println("cover_lp["+i+"] = "+cover_lp[i]);
		}
		
		for (int i =0;i<Num_user;i++){
			double max_val =0;
			int max_index = 0;
			for (int j=0;j<V_num_AP;j++){
				if(cover_lp[j*Num_user+i]>0){
					int real_index = get_real_index(j,Num_Ch);
					if(Double.compare(Weight[real_index][i]*cover_lp[j*Num_user+i], max_val)==1){ // maximum case
						max_val = Weight[real_index][i]*cover_lp[j*Num_user+i];
						max_index = j;
					}
				}
			}
			int real_index = get_real_index(max_index,Num_Ch);
			fix_sol[real_index][i] =1;
			operating_channel[real_index] =get_channel_index (max_index,Num_Ch);
		}
		return fix_sol;
			
	}
	
	// find an optimal solution of LP using CPLEX
	// return the solution from CPLEX
	public double [] set_cover_rounding(int [] V_on){
		
		double [] LP_sol = new double [Num_user*Num_AP*Num_Ch+Num_AP*Num_Ch]; // solution of LP
		try {
			int d_max = Cal_d_max(Num_AP,Num_user,Weight);
			int[][] V_I_AP = virtual_interference_graph(Overlap_G,Num_Ch,Num_AP);
			
			IloCplex cplex = new IloCplex();
			
			// association variable & on_off variable
			// z11 z21 .. zU1 z12 z22 .. zU2 .. z1N z2N .. zUN y11 y12 ..y1c y21 y22 .. y2c .. yN1 yN2 .. yNc
			IloNumVar [] z = cplex.numVarArray(Num_user*Num_AP*Num_Ch + Num_AP*Num_Ch, 0.0, 1);
		
			
			// set the coefficient for z
			double [] W = new double[Num_user*Num_AP*Num_Ch+Num_AP*Num_Ch];		
			
			for (int j =0;j<Num_AP;j++){
				for (int c = 0; c<Num_Ch;c++){
					W[Num_user*Num_AP*Num_Ch + j*Num_Ch+c] = -Cost_AP[j];
					for(int i = 0;i<Num_user;i++){
						if(V_on[j*Num_Ch+c]==1){
							W[j*Num_Ch*Num_user+Num_user*c+i]= Weight[j][i];
						}
						System.out.println(W[j*Num_Ch*Num_user+Num_user*c+i]);
					}
				}
			}
			
			// objective fn
			IloLinearNumExpr objective = cplex.scalProd(z, W);
			
			// max obj
			cplex.addMaximize(objective);
			
			// constraints 
			// single association & cover all users
			if(d_max ==0){
				return LP_sol;
			}
			double epsilon = (double) 1/(2*d_max);
			for (int i =0;i<Num_user; i++){
				double [] coeff = new double[Num_user*Num_AP*Num_Ch+ Num_AP*Num_Ch];
				int temp_flag =0;
				for (int j =0;j<Num_AP*Num_Ch;j++){
					if(W[i+j*Num_user]>0){
						coeff[i+j*Num_user] =1;
						temp_flag =1;
					}
				}
				if(temp_flag==1){
				// sum_{j} zij <1(lambda)
				cplex.addLe(cplex.scalProd(z, coeff),1);
				// sum_{j} zij > epsilon
				System.out.println("epsilon : " + epsilon + "  d_max : " + d_max);
				cplex.addGe(cplex.scalProd(z, coeff),epsilon);
				}
			}
			
			//configuration 
			for (int j=0;j<Num_AP*Num_Ch;j++){
				double [] coeff2 = new double[Num_user*Num_AP*Num_Ch+ Num_AP*Num_Ch];
				
				for (int i=0; i<Num_user;i++){
					coeff2[j*Num_user +i]=1;
				}
				coeff2[Num_user*Num_AP*Num_Ch+j] =-1;
				
				//sum_{i} zij < sum_{c} yjc
				cplex.addLe(cplex.scalProd(z, coeff2),0);
				
			}
			//interference
			for (int j=0;j<Num_AP*Num_Ch;j++){
				for (int i =0; i<Num_AP*Num_Ch;i++){
					if(V_I_AP[j][i]==1){
						double [] coeff3 = new double[Num_user*Num_AP*Num_Ch+ Num_AP*Num_Ch];
						coeff3[Num_user*Num_AP*Num_Ch+i] =1;
						coeff3[Num_user*Num_AP*Num_Ch+j] =1;
						cplex.addLe(cplex.scalProd(z, coeff3),1);
					}						
				}
			}
			 
			// solve
			if(cplex.solve()){
				System.out.println("model solved");
				System.out.println("obj = "+cplex.getObjValue());
				LP_sol = cplex.getValues(z);
			}else{
				System.out.println("model not solved");
			}
			
		} catch (IloException exc) {
			exc.printStackTrace();
		}
		
		return LP_sol;
	}
		
		
}
