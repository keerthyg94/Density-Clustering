import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.lang.Math;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;
import java.util.TreeMap;

public class density {

	int g=2;
	static int row1=0;
	static int column1=0;
	static int mGenes=0;
	static int groundtruth[];
	static int foundtruth[];
	static double data_matrix[][] = new double[600][16];
	static double distance_matrix [] []; 
	static double new_distance_matrix[][];
	static boolean visited[];
	static ArrayList<Integer> Neighborpts = new ArrayList<Integer>(); 
	static Map<Integer,Integer> Clustermap = new TreeMap<Integer,Integer>();
    

	/*
	 * DBSCAN(D, eps, MinPts)
       C = 0
       for each unvisited point P in dataset D
          mark P as visited
          NeighborPts = regionQuery(P, eps)
          if sizeof(NeighborPts) < MinPts
          mark P as NOISE
           else
          C = next cluster
          expandCluster(P, NeighborPts, C, eps, MinPts)
	 */
	static void DBSCAN(double eps,int minpts)
	{	

		int cluster = 0;

		for (int p = 0; p < mGenes;p++)
		{
			//System.out.println("Gene"+p+" "+visited[p-1]);
			if(!visited[p])
			{
				visited [p] = true;

				Neighborpts = regionQuery(p,eps);

				if (Neighborpts.size() < minpts)
				{
					//System.out.println(p+"is noise");
					Clustermap.put(p,-1);
				}
				else 
				{
					cluster ++;
					//System.out.println((p+1)+"gene in cluster"+cluster);
					expandCluster(p,Neighborpts,cluster,eps,minpts);
				}
			}
		}
		
		System.out.println("Number of clusters "+cluster);
		System.out.println();

	}


	/*
	 * regionQuery(P, eps)
	 * return all points within P's eps-neighborhood
	 */

	static ArrayList<Integer> regionQuery(int p, double eps)
	{
		ArrayList<Integer> neigh = new ArrayList<Integer>();
		int gindex = p;

		for (int i=0;i<mGenes;i++)
		{
			if (distance_matrix[gindex][i] <= eps)
			{
				if(i!=gindex)
				{
					neigh.add(i);
				}
			}
		}
		//System.out.println(neigh);
		return neigh;
	}

	/*
	 * expandCluster(P, NeighborPts, C, eps, MinPts)
       add P to cluster C
   		for each point P' in NeighborPts 
      	if P' is not visited
         mark P' as visited
         NeighborPts' = regionQuery(P', eps)
         if sizeof(NeighborPts') >= MinPts
            NeighborPts = NeighborPts joined with NeighborPts'
      	if P' is not yet member of any cluster
         add P' to cluster C
	 */

	static void expandCluster(int p,ArrayList<Integer> neighbor,int cluster,double eps,int minpts)
	{
		Clustermap.put(p,cluster);
		ArrayList<Integer> newneigh = new ArrayList<Integer>();

		//System.out.println(neighbor.size());

		for (int i=0;i<neighbor.size();i++)
		{
			//System.out.println("size"+neighbor.size());
			int p1 = neighbor.get(i);
			if( visited[p1] != true)
			{
				visited [p1] = true;
				newneigh = regionQuery(p1,eps);
				if(newneigh.size()>=minpts)
				{
					expandCluster(p1,newneigh,cluster,eps,minpts);
					Clustermap.put(p1, cluster);
				
				}
				else
				{
					Clustermap.put(p1, cluster);
					//
				}

			}
		}

	}

	public static void main(String[] args)
	{

		double eps= 0;
		int minpts = 0;

		
		String filename=null;
		System.out.println("Density based clustering");
		System.out.println();
		System.out.println("Enter the file");
		try{
			BufferedReader reader=new BufferedReader(new InputStreamReader(System.in));
			filename = reader.readLine();
			
			if (filename.equals("iyer.txt"))
			{
				eps=1.2;
				minpts=1;
		         
			}
			else if (filename.equals("cho.txt"))
			{
				eps= 1.24;
				minpts = 7;
			}


			BufferedReader readfile = new BufferedReader(new FileReader(filename));
			while(readfile.readLine()!=null)	
			{
				mGenes++;
			}
			System.out.println("Total no of genes ="+mGenes);
			
			@SuppressWarnings("resource")
			BufferedReader readfile1 = new BufferedReader(new FileReader(filename));
			String read;
			read=readfile1.readLine();
			String tokens[]=read.split("\\t");
			int columns =tokens.length-2;
			System.out.println();
			System.out.println(" Number of columns = "+columns);
			groundtruth= new int[mGenes];
			BufferedReader readfile2 = new BufferedReader(new FileReader(filename));

			String read1;

			while ((read1 = readfile2.readLine()) != null) 
			{
				String tokens1[]=read1.split("\\t");
				column1=0;
				for(int i=0;i<=1;i++)
				{
					String geneid=tokens1[0];
					String truth=tokens1[1];
					int gid=Integer.parseInt(geneid);
					int gtruth=Integer.parseInt(truth);
					groundtruth[row1]= gtruth;
				}
				for(int i=2;i<tokens1.length;i++)
				{
					double value=Double.parseDouble(tokens1[i]);
					data_matrix[row1][column1]=value;
					column1++;

				}
				row1++;
			}

			visited = new boolean[mGenes];
			for(int i=0;i<mGenes;i++)
			{
				visited[i]=false;
			}
			
			distance_matrix= new double[mGenes][mGenes];
			new_distance_matrix = new double[mGenes][mGenes];
			foundtruth = new int[mGenes];
			//System.out.println();
			for(int i=0;i<row1;i++)
			{
				for(int j=0;j<row1;j++)
				{
					if(i==j)
					{
						distance_matrix[i][j]=0.0;
						new_distance_matrix[i][j]=0.0;
					}
					else 
					{	double sum=0;
					double temp=0;
					for(int k=0;k<column1;k++)
					{  
						temp=data_matrix[i][k]-data_matrix[j][k];
						double sq_temp=temp*temp;
						sum=sum+sq_temp;
					}
					double root_temp=Math.sqrt(sum);
					distance_matrix[i][j]=root_temp;
					new_distance_matrix[i][j]=root_temp;

					}
				}

			}
			
		}
	
		catch(IOException e)
		{
			System.out.println();
		}

		DBSCAN(eps,minpts);
		int count=0;
		for(Map.Entry<Integer, Integer> entry: Clustermap.entrySet())
		{
			//if(entry.getValue()!=71)
			foundtruth[entry.getKey()]=entry.getValue();
			//System.out.println((entry.getKey()+1) +"\t" + (entry.getValue()));
			System.out.println(" Gene Id = "+(entry.getKey()+1) + " Cluster number = "+entry.getValue());
			count++;
		}
		//System.out.println(Clustermap);
		System.out.println("Size of the cluster map = " +count); 
		
		//Incident_matrix calculation
				int incident_matrixP[][] = new int[mGenes][mGenes];
				int incident_matrixC[][] = new int[mGenes][mGenes];
				int ss=0,dd=0,sd=0,ds=0;
				for(int i=0;i<mGenes;i++)
				{
					for(int j=0; j<mGenes; j++)
					{
						if(groundtruth[i]==groundtruth[j])
						{
							incident_matrixP[i][j] =1;
						}
						
						else
						{
							incident_matrixP[i][j]=0;
						}
						
						if(foundtruth[i]==foundtruth[j])
						{
							incident_matrixC[i][j]=1;
						}
						else
						{
							incident_matrixC[i][j]=0;
						}
						
					}
				}//end of for
				
			for(int i=0;i<mGenes;i++)	
			{
				for(int j=0;j<mGenes;j++)
				{
					if(incident_matrixP[i][j]==incident_matrixC[i][j])
					{
						if(incident_matrixP[i][j]==1)
						{
							ss++;
						}
						else
						{
							dd++;
						}
					}
					else
					{
						if(incident_matrixP[i][j]==1)
						{
							ds++;
						}
						else
						{
							sd++;
						}
					}
				}
			}
				
			System.out.println("SS = "+ss);
			System.out.println("DD = "+dd);
			System.out.println("SD = "+sd);
			System.out.println("DS = "+ds);
			
			double rand_index;
			double jaccard_coef;
			double ss1=(double)ss;
			double dd1 = (double)dd;
			double sd1= (double)sd;
			double ds1 = (double)ds;
			
			rand_index =(ss1+dd1)/(ss1+dd1+sd1+ds1);
			jaccard_coef=(ss1)/(ss1+sd1+ds1);
			
			System.out.println(" Rand index = "+rand_index);
			System.out.println(" Jacard Coefficient = "+jaccard_coef);
			
			//internal_index calculation
			double cincident_matrix[][] = new double[mGenes][mGenes];
			for(int i=0;i<mGenes;i++)
			{
				for(int j=0; j<mGenes;j++)
				{
					cincident_matrix[i][j] =(double)incident_matrixC[i][j];
				}
			}
			
			double dmean = mean(new_distance_matrix,mGenes);
			double cmean = mean(cincident_matrix,mGenes);
			double nvar=numerical_variance(new_distance_matrix, cincident_matrix, mGenes);
			//System.out.println(" Nvar = "+nvar);
			double dvar=Math.sqrt(variance(new_distance_matrix,mGenes));
			//System.out.println(" Dvar = "+dvar);
			double cvar = Math.sqrt(variance(cincident_matrix, mGenes));
			//System.out.println(" Cvar = "+cvar);
			double correlation = Math.abs(nvar/(dvar*cvar));
			System.out.println(" Correlation of incident matrix and distance matrix = "+correlation);
	} 
	
	public static double mean(double arr[][], int num)
	{
		double sum=0,mean=0;
		for(int i=0;i<num;i++)
		{
			for( int j=0;j<num;j++)
			{
				sum=sum+arr[i][j];
			}
		}
		mean = sum/(num*num);
			
		return mean;
	}
	
	public static double variance(double arr[][], int num)
	{
		double sum=0,mean=0, var=0,sum_mean=0;
		for(int i=0;i<num;i++)
		{
			for( int j=0;j<num;j++)
			{
				sum=sum+arr[i][j];
			}
		}
		mean = sum/(num*num);
		for(int i=0;i<num;i++)
		{
			for( int j=0;j<num;j++)
			{
				sum_mean= sum_mean + ((arr[i][j]-mean)*(arr[i][j]-mean));
			}
		}
		return sum_mean;
	}
	
	public static double numerical_variance(double arr[][],double arr1[][],int num)
	{
		double sum=0,sum1=0,mean=0,mean1=0,sum_mean=0;
		for(int i=0;i<num;i++)
		{
			for( int j=0;j<num;j++)
			{
				sum=sum+arr[i][j];
				sum1=sum1+arr1[i][j];
			}
		}
		mean = sum/(num*num);
		mean1=sum1/(num*num);
		for(int i=0;i<num;i++)
		{
			for( int j=0;j<num;j++)
			{
				sum_mean= sum_mean + ((arr[i][j]-mean)*(arr1[i][j]-mean1));
			}
		}
		return sum_mean;
	}
	
	public static int num_elements(int[] arr)
	{
		Set<Integer> newset = new HashSet<Integer>();
		for(int element: arr)
		{
			newset.add(element);
		}
		if(newset.contains(-1))
		{
			return newset.size()-1;
		}
		else
		{
			return newset.size();
		}
		
	}
}
