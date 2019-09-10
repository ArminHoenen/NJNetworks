package mainpack;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

public class NJNetworks {
	
	public static ArrayList<String> njtrees;
	//the NJ network is hierarchically built bottom-up; 
	//first i generate all nj trees
	public static String[] initialTaxa;
	public static String mode;
	public static int digit;
	
	public static void main (String [] args) throws Exception{
		
		//TODO
		//a) identify mildest rounding factor b) increment rounding c) undulate mode halfup/down d)properly output rounding params
				
		mode="half-up";
		digit=0;
		
		double [][] mat=getMat1(new File("dmatred.csv"));//args[0]
		String [] taxa = getTaxa1(new File("dmatred.csv"));
		//for more than 14 taxa (with mathalf 88 cells)   100 below must be incremented
		initialTaxa=taxa;
		
		
		boolean stop=false;
		boolean printed=false;
		
		while (!stop) {
			
			ArrayList<String[]>restrict= new ArrayList<String[]>();
			njtrees=new ArrayList<String>();
			ArrayList<String[]> njtree=new ArrayList<String[]>();
			
			nextNJStep(mat,njtree,restrict,taxa);			
			
			if (njtrees.size()>1) {
				
				if (digit==0) {
					System.out.println(njtrees.size()+" naturally occurring nj trees:");
				}else {
					System.out.println("rounding to next integer in mode "+mode);
					System.out.println(njtrees.size()+" nj trees:");
				}
				for (String at:njtrees) {
					System.out.println(at);
				}
				stop=true;
			}
			else {
				
				digit+=1;
				if (!printed) {
					System.out.println("Only one natural NJ tree with matrix:\n"+njtrees.get(0)+"\nNJN initiates rounding...");
					printed=true;
				}			
			}

		}//while ! stop
		
		
		
		
		
		
		
	}//main
	
public static double[][] getMat1(File in)throws Exception{
		
		String line="";
		BufferedReader br =  new BufferedReader(
				 new InputStreamReader(
						 new FileInputStream(in), "UTF-8"));

		ArrayList<String []> dat= new ArrayList<String[]>();
		while ((line=br.readLine())!=null){
			
			dat.add(line.split("\t"));
			
		}
		br.close();
		
		double [][] res= new double[dat.size()-1][dat.size()-1];
		
		int rowind=0;
		for (int a=1;a<dat.size();a++) {
			
			String [] act=dat.get(a);
			for (int b=1;b<act.length;b++) {
				
				res[rowind][(b-1)]=Double.parseDouble(act[b]);
				
			}
			rowind+=1;
			
		}

		return res;
		}//getMat
public static String[] getTaxa1(File in)throws Exception{
	
	String line="";
	BufferedReader br =  new BufferedReader(
			 new InputStreamReader(
					 new FileInputStream(in), "UTF-8"));

	String[] res=new String[] {};
	while ((line=br.readLine())!=null){
		res=line.split("\t");
		break;
		
	}
	br.close();

	return res;
	}//getTaxa
	
	
	public static void nextNJStep (double[][] mat, ArrayList<String[]> tree, ArrayList<String[]> restrict, String[] taxa) throws Exception{
		
		
		ArrayList<String[]> treeOne=(ArrayList<String[]>)tree.clone();
		
		//System.out.println("*********************STEP");
		if (taxa.length>1) {
			
			//compute mediate matrix
			double [][] medmat=new double[taxa.length][taxa.length];
			
			double [] avev=new double[taxa.length];
			
			double res=0;
			for (int a=0;a<taxa.length;a++) {
				for (int b=0;b<taxa.length;b++)res+=mat[a][b];
				if (taxa.length>2)res/=(double)(taxa.length-2);//if is for last step, when the divergence concerns only the 2 remaining taxa
				//if (res==Double.POSITIVE_INFINITY)System.out.println(taxa.length); >> 2
				avev[a]=res;
			}
			for (int a=0;a<taxa.length;a++)medmat[a][a]=0;
			for (int a=0;a<taxa.length;a++) {
				for (int b=0;b<taxa.length;b++) {
					if (a<b) {
						double val=mat[a][b]-(avev[a]+avev[b]);
						//System.out.println(val+" "+Math.round(val));
						//System.out.println(val+" "+roundAvoid((mat[a][b]-(avev[a]+avev[b])),2));
						if (digit==0) {
							medmat[a][b]=val;//val
							medmat[b][a]=val;//val
						}
						else {
						medmat[a][b]=Math.round(val);//val
						medmat[b][a]=Math.round(val);//val
						}
						//medmat[a][b]=roundAvoid(val,1);//val
						//medmat[b][a]=roundAvoid(val,1);//val
					}
				}
			}
			//System.out.println("MedMat:");
			//printMat(medmat,taxa.length);
			
			//look for smallest val > limit
			double min=1000.0;//101
			//find min
			for (int a=0;a<taxa.length;a++) {
				for (int b=0;b<taxa.length;b++) {
					if (a<b) {
						if (medmat[a][b]<min) {
							min=medmat[a][b];
						}
						
					}
				}
		      }//find new min
			//System.out.println(min);
			
			//find next 
			for (int a=0;a<taxa.length;a++) {
				for (int b=0;b<taxa.length;b++) {
					if (a<b) {
						if (medmat[a][b]==min) {
						//check if (a,b) in restriction set, otherwise take, add to restriction set and shoot new path
						//System.out.println(taxa[a]+taxa[b]);
						
							//System.out.println("*******not is in");
							String[] act=new String[] {taxa[a],taxa[b]};
											
								//update matrix step 2 into mat
								treeOne.add(act);
								
								double [][] newMat=new double[taxa.length+1][taxa.length+1];
								for (int c=0;c<taxa.length+1;c++) newMat[c][c]=0;
								for (int c=0;c<taxa.length+1;c++) {
									
									for (int d=0;d<taxa.length+1;d++) {
										
										if (c<d) {
											
											if (d!=taxa.length) {
												newMat[c][d]=mat[c][d];
												newMat[d][c]=mat[c][d];
											}
											else {
												//mema
												double newVal=0;
												newVal=(mat[c][a]+mat[c][b])-mat[a][b];
												newVal/=2;
												newMat[c][d]=newVal;
												newMat[d][c]=newVal;										
											}
											
										}
									}
									
								}
								//System.out.println("NewMat:");
								//printMat(newMat,taxa.length+1);
								double[][] newnewMat=new double[taxa.length-1][taxa.length-1];
																
								int index1=0;
								int index2=0;
								
								for (int x=0;x<taxa.length+1;x++) {
									index2=0;
									
									if (x!=a&&x!=b) {
										
										for (int y=0;y<taxa.length+1;y++) {
											
											if (y!=a&&y!=b) {
											newnewMat[index1][index2]=newMat[x][y];
											newnewMat[index2][index1]=newMat[x][y];
											index2+=1;
											}
										}
										index1+=1;
									}
									
								}
								//System.out.println("NewNewMat:");
								//printMat(newnewMat,taxa.length-1);
								
								//new taxa
								ArrayList<String> tp=new ArrayList<String>();
								for (int x=0;x<taxa.length;x++) {
									if (x!=a&&x!=b)tp.add(taxa[x]);
								}
								String newTax=taxa[a]+"-"+taxa[b];
								tp.add(newTax);
								String[] tpp=new String[taxa.length-1];
								for (int z=0;z<tp.size();z++)tpp[z]=tp.get(z);
										
								restrict = new ArrayList<String[]>();
								nextNJStep (newnewMat, treeOne, restrict, tpp);
								treeOne=(ArrayList<String[]>)tree.clone();
													
								
						}
					}
				}
			}
			
			
			
			//check restriction set
			
			//if restriction set exhausted, nextNJStep level up
			
			
			
		}//nj tree not complete
		else {
			//transform data structure into other format:
			//		
			
			//System.out.print("**NEW tree: ");
			//for (String [] ata:tree) {
			//	System.out.print("("+ata[0]+","+ata[1]+")");
			//}
			//System.out.print("\n");
			String tree_reformatted=getReformatted(tree);
			//System.out.print("\n");
			if (!njtrees.contains(tree_reformatted))njtrees.add(tree_reformatted);
		}
		
		
	}
	
	public static String getReformatted (ArrayList<String[]> in_tree) throws Exception{
		
		//order all alphabetically and for level (number of nodes connected) > per tree only one order
		//all that have 2, then 3, ... order inner, order outer
		String res="";
		
		for (int an = 2;an<initialTaxa.length+1;an++) {
			
			ArrayList<String []> level= new ArrayList<String[]>();
			
			for (String [] act:in_tree) {
				int num=0;
				String[] els1=act[0].split("-");
				String[] els2=act[1].split("-");
				num+=els1.length;
				num+=els2.length;
				
				if (num==an)level.add(act);
				
			}
			
			if (level.size()>0) {
				
				ArrayList<String []> l_oi= new ArrayList<String[]>();//level_ordered_inner
				//order inner
				for (String [] act:level) {
					String [] oi=new String[2];
					String el1=getAlphabetical(act[0]);
					String el2=getAlphabetical(act[1]);
					
					if (el1.charAt(0)<el2.charAt(0)) {
						oi[0]=el1;
						oi[1]=el2;
					}
					else {
						oi[0]=el2;
						oi[1]=el1;
					}
					l_oi.add(oi);
				}
				
				//l_oi now has String arrays with ordered pairs
				//order outer
				String firsts ="";
				for (int a=0;a<l_oi.size()-1;a++) {
					String aoi=l_oi.get(a)[0];
					firsts+=aoi.charAt(0)+"-";
				}
				firsts+=l_oi.get(l_oi.size()-1)[0].charAt(0);
				String f_o=getAlphabetical(firsts);
				
				String [] els=f_o.split("-");
				
				for (int b=0;b<els.length;b++) {
					for (int a=0;a<l_oi.size();a++) {
						if (l_oi.get(a)[0].startsWith(els[b])) {
							res+="("+l_oi.get(a)[0]+","+l_oi.get(a)[1]+")";
						}
					}
				}
				
			}//if level has input
			
			
		}//for levels
		
		//System.out.println("RES ordered: "+res);
		
		return res;
	}
	
	public static String getAlphabetical(String in) throws Exception{
		
		//System.out.println("**IN**"+in);
		
		String res="";
		String [] els=in.split("-");
		String newin ="";
		for (String el:els)newin+=el;
		
		char[] charArray = newin.toCharArray();
		Arrays.sort(charArray);
		String sortedString = new String(charArray);
		//System.out.println(sortedString);  
		
		for (int a=0;a<sortedString.length()-1;a++) {
			res+=sortedString.charAt(a)+"-";
		}
		res+=sortedString.charAt(sortedString.length()-1);
		
		//System.out.println("**OUT"+res);
		
		return res;
	}
	
	
	public static boolean hasAll(double[][] mat, double limit, ArrayList<String[]> res, String[] taxa) throws Exception{
		
		for (int a=0;a<taxa.length;a++) {
			for (int b=0;b<taxa.length;b++) {
				
				if (a<b) {
				  if (mat [a][b]==limit) {
					  if (!isIn(taxa[a],taxa[b],res))return false;
				  }
				}
					
			}
		}
		
		return true;
	}
	public static boolean isIn(String a, String b, ArrayList<String[]> res) throws Exception{
		
		for (String [] ac:res) {
			
			if ((ac[0].equals(a)&&ac[1].equals(b))||(ac[0].equals(b)&&ac[1].equals(a)))return true;
			
		}
		return false;
		
	}
	
	
	public static void printMat(double [][] mat, int taxa)throws Exception{
		for (int a=0;a<taxa;a++) System.out.print("\""+a+"\"\t");
	       System.out.print("\n");
			for (int a=0;a<taxa;a++) {
				System.out.print("\""+a+"\"\t");
				for (int b=0;b<taxa;b++) {
					System.out.print(mat[a][b]+"\t");
				}
				System.out.print("\n");
			}
	}
	public static double roundAvoid(double value, int places) {
	    double scale = Math.pow(10, places);
	    return Math.round(value * scale) / scale;
	}

}
