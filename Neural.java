import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Random;
import java.util.Scanner;
///////////////////////////////////////////////////////////////////////////////
//ALL STUDENTS COMPLETE THESE SECTIONS
//Title:            (Neural)
//Files:            (Neural.java)
//Semester:         (CS540) Fall 2018
//
//Author:           (Jermiah Ansera)
//Email:            (ansera@wisc.edu)
//CS Login:         (jermiah)
//Lecturer's Name:  (Jerry Zhu)
////////////////////PAIR PROGRAMMERS COMPLETE THIS SECTION ////////////////////
class setItem {

    public double X1;
    public double X2;
    public double Y;

    public setItem(double x1, double x2, double y) {
        X1 = x1;
        X2 = x2;
        this.Y = y;
    }

}
public class Neural {
	 private static double[] WEIGHTS;
	 private static double X1;
	 private static double X2;
	 private static double Y;
	 private static double eta;
	 private static final String EVALUATION = "hw2_midterm_A_eval.txt";
	 private static final String TESTING = "hw2_midterm_A_test.txt";
	 private static final String TRAINING = "hw2_midterm_A_train.txt";   
	 
	 private static double[] FinduAvAuBvBuCvC() {
	        double uA = WEIGHTS[0] + ( WEIGHTS[1] * X1 ) + (WEIGHTS[2] * X2);
	        double vA = Math.max(uA, 0); //ReLu
	        
	        double uB = WEIGHTS[3] + (WEIGHTS[4] * X1) + (WEIGHTS[5] * X2);
	        double vB = Math.max(uB, 0); //ReLu
	        
	        double uC = WEIGHTS[6] + (WEIGHTS[7] * vA) + (WEIGHTS[8] * vB);
	        
	        double vC = 1.0 / (1 + Math.exp(0 - uC)); //SigmoidC
	        
	        return new double[] {uA, vA, uB, vB, uC, vC};
	    }
	 private static double[] DerivativesOfOutputlayer(double vC) {
		 	double E = (Math.pow(vC - Y, 2)) / 2;
		 	double dEdvC = vC - Y;
		 	double dEduC = dEdvC * (vC * (1-vC));
	        return new double[] {E, dEdvC, dEduC};
	    }
	 private static double[] DerivativesOfHiddenLayer(double vC, double uA, double uB) {
		 double[] output = DerivativesOfOutputlayer(vC);
		 	double dEduK = output[2];
		 
		 	double dEdvA = WEIGHTS[7] * dEduK;
		 	double dEduA; 
		 	
		 	if (uA >= 0) {
		 		dEduA = dEdvA; 
		 	} else {
		 		dEduA = 0;
		 	}
		 	
		 	double dEdvB = WEIGHTS[8] * dEduK;
		 	double dEduB; 
		 	
		 	if (uB >= 0) {
		 		dEduB = dEdvB;
		 	} else {
		 		dEduB = 0;
		 	}
		 	
	        return new double[] {dEdvA, dEduA, dEdvB, dEduB};
	    }
	 private static double[] DerivativesOfWeights(double[] HiddenLayer, double[] output, 
			 double vA, double vB) {
		 double[] dEdW = new double[WEIGHTS.length];
		 //weights for node A
		 dEdW[0] = 1 * HiddenLayer[1];
		 dEdW[1] = X1 * HiddenLayer[1];
		 dEdW[2] = X2 * HiddenLayer[1];
		 //weights for node B
		 dEdW[3] = 1 * HiddenLayer[3];
		 dEdW[4] = X1 * HiddenLayer[3];
		 dEdW[5] = X2 * HiddenLayer[3];
		 
		 dEdW[6] = 1 * output[2];
		 dEdW[7] = vA * output[2];
		 dEdW[8] = vB * output[2];
		 
		 //how do i get rid of -0.0???
		 for(int i = 0; i < dEdW.length; i++) {
			 dEdW[i] = Math.round(dEdW[i] * 10000000.0) / 10000000.0;
		 }
		 		 
		 return dEdW;
	 }
	 private static void updateWEIGHTSsgd(double[] dEdW) {
		 for (int i = 0; i < WEIGHTS.length; i++) {
	            WEIGHTS[i] = WEIGHTS[i] - (eta * dEdW[i]);
	        }
	 }
	 public static ArrayList<setItem> FileReadSetItems(String filename) {
	        ArrayList<setItem> setItems = new ArrayList<>();
	        File file = new File(filename);
	        Scanner sc;
			try {
				sc = new Scanner(file);
				while (sc.hasNextLine()) {
	                String line = sc.nextLine();
	                String[] tokens = line.split("\\s+");
	                setItems.add(new setItem(Double.parseDouble(tokens[0]),
	                		Double.parseDouble(tokens[1]), Double.parseDouble(tokens[2])
	                		));
	            }
			} catch (FileNotFoundException e) {
				System.out.println("File " + filename + " is not found.");
	            System.exit(1);
			}
	       	       
	        return setItems;
	    }
	 private static double setError(ArrayList<setItem> SetItems) {
		 double Error = 0.0;
		 for (int i = 0; i < SetItems.size(); i++) {			 
			X1 = SetItems.get(i).X1;
			X2 = SetItems.get(i).X2;				
			Error = Error + Math.pow((FinduAvAuBvBuCvC()[5] - SetItems.get(i).Y), 2);
		 }
		 
		 Error = Error * 0.5;
		 
		 return Error;
	 }
	 
	 
	public static void main(String[] args) {
		
		try {
			int flag = Integer.parseInt(args[0]); //flag number
			
			
			if(flag == 100) {
				WEIGHTS = new double[9];
				for (int i = 0; i < WEIGHTS.length; i++) {
	                WEIGHTS[i] = Double.parseDouble(args[i + 1]);
	            }
				
				X1 = Double.parseDouble(args[10]);
				X2 = Double.parseDouble(args[11]);
				
				double[] uAvAuBvBuCvC = new double[6];
				uAvAuBvBuCvC = FinduAvAuBvBuCvC();
				System.out.printf( "%.5f ", uAvAuBvBuCvC[0]);	
				System.out.printf( "%.5f ", uAvAuBvBuCvC[1]);	
				System.out.printf( "%.5f ", uAvAuBvBuCvC[2]);	
				System.out.printf( "%.5f ", uAvAuBvBuCvC[3]);	
				System.out.printf( "%.5f ", uAvAuBvBuCvC[4]);	
				System.out.printf( "%.5f", uAvAuBvBuCvC[5]);	
			}
			else if (flag == 200) {
				WEIGHTS = new double[9];
				for (int i = 0; i < WEIGHTS.length; i++) {
	                WEIGHTS[i] = Double.parseDouble(args[i + 1]);
	            }
				
				X1 = Double.parseDouble(args[10]);
				X2 = Double.parseDouble(args[11]);
				double[] uAvAuBvBuCvC = new double[6];
				uAvAuBvBuCvC = FinduAvAuBvBuCvC();
				Y = Double.parseDouble(args[12]);
				
				double[] output = new double[3];
				output = DerivativesOfOutputlayer(uAvAuBvBuCvC[5]);

				System.out.printf( "%.5f ", output[0]);	
				System.out.printf( "%.5f ", output[1]);	
				System.out.printf( "%.5f", output[2]);				
			}
			else if (flag == 300) {
				WEIGHTS = new double[9];
				for (int i = 0; i < WEIGHTS.length; i++) {
	                WEIGHTS[i] = Double.parseDouble(args[i + 1]);
	            }
				
				X1 = Double.parseDouble(args[10]);
				X2 = Double.parseDouble(args[11]);
				double[] uAvAuBvBuCvC = new double[6];
				uAvAuBvBuCvC = FinduAvAuBvBuCvC();
				Y = Double.parseDouble(args[12]);
				
				double[] output = new double[3];
				output = DerivativesOfHiddenLayer(uAvAuBvBuCvC[5], uAvAuBvBuCvC[0], uAvAuBvBuCvC[2]);

				System.out.printf( "%.5f ", output[0]);	
				System.out.printf( "%.5f ", output[1]);	
				System.out.printf( "%.5f ", output[2]);	
				System.out.printf( "%.5f", output[3]);	
			} 
			else if (flag == 400) {
				WEIGHTS = new double[9];
				for (int i = 0; i < WEIGHTS.length; i++) {
	                WEIGHTS[i] = Double.parseDouble(args[i + 1]);
	            }
				
				X1 = Double.parseDouble(args[10]);
				X2 = Double.parseDouble(args[11]);
				
				double[] uAvAuBvBuCvC = new double[6];
				uAvAuBvBuCvC = FinduAvAuBvBuCvC();
				
				Y = Double.parseDouble(args[12]);
				
				double[] HiddenLayer = new double[3];
				HiddenLayer = DerivativesOfHiddenLayer(uAvAuBvBuCvC[5], uAvAuBvBuCvC[0], uAvAuBvBuCvC[2]);
				
				double[] output = new double[3];
				output = DerivativesOfOutputlayer(uAvAuBvBuCvC[5]);
				
				double[] derivativesOfWeights = DerivativesOfWeights(HiddenLayer, output, uAvAuBvBuCvC[1], uAvAuBvBuCvC[3]);
				
				System.out.printf( "%.5f ", derivativesOfWeights[0]);	
				System.out.printf( "%.5f ", derivativesOfWeights[1]);	
				System.out.printf( "%.5f ", derivativesOfWeights[2]);	
				System.out.printf( "%.5f ", derivativesOfWeights[3]);	
				System.out.printf( "%.5f ", derivativesOfWeights[4]);	
				System.out.printf( "%.5f ", derivativesOfWeights[5]);	
				System.out.printf( "%.5f ", derivativesOfWeights[6]);	
				System.out.printf( "%.5f ", derivativesOfWeights[7]);
				System.out.printf( "%.5f", derivativesOfWeights[8]);	
				
			} 
			else if (flag == 500) {
				WEIGHTS = new double[9];
				for (int i = 0; i < WEIGHTS.length; i++) {
	                WEIGHTS[i] = Double.parseDouble(args[i + 1]);
	            }
				
				X1 = Double.parseDouble(args[10]);
				X2 = Double.parseDouble(args[11]);
				double[] uAvAuBvBuCvC = new double[6];
				uAvAuBvBuCvC = FinduAvAuBvBuCvC();
				
				Y = Double.parseDouble(args[12]);
				eta = Double.parseDouble(args[13]);
						
				double[] output = new double[3];
				output = DerivativesOfOutputlayer(uAvAuBvBuCvC[5]);
				
				for (int i = 0; i < WEIGHTS.length - 1; i++) {
					System.out.printf( "%.5f ", WEIGHTS[i]);
	            }
				System.out.printf( "%.5f", WEIGHTS[WEIGHTS.length - 1]);
				System.out.printf( "\n%.5f\n", output[0]);	
				
				double[] HiddenLayer = new double[3];
				HiddenLayer = DerivativesOfHiddenLayer(uAvAuBvBuCvC[5], uAvAuBvBuCvC[0], 
						uAvAuBvBuCvC[2]);
								
				double[] derivativesOfWeights = DerivativesOfWeights(HiddenLayer, output, 
						uAvAuBvBuCvC[1], uAvAuBvBuCvC[3]);
				
				updateWEIGHTSsgd(derivativesOfWeights);
				for (int i = 0; i < WEIGHTS.length - 1; i++) {
					System.out.printf( "%.5f ", WEIGHTS[i]);
	            }
				System.out.printf( "%.5f", WEIGHTS[WEIGHTS.length - 1]);
				
				uAvAuBvBuCvC = FinduAvAuBvBuCvC();
				output = DerivativesOfOutputlayer(uAvAuBvBuCvC[5]);
				System.out.printf( "\n%.5f\n", output[0]);
			} 
			else if (flag == 600) {
				eta = Double.parseDouble(args[10]);
				
				WEIGHTS = new double[9];
				for (int i = 0; i < WEIGHTS.length; i++) {
	                WEIGHTS[i] = Double.parseDouble(args[i + 1]);
	            }
				ArrayList<setItem> training = FileReadSetItems(TRAINING);
				ArrayList<setItem> evaluation = FileReadSetItems(EVALUATION);
				for (int i = 0; i < training.size(); i++) {
					X1 = training.get(i).X1;
					X2 = training.get(i).X2;
					Y = training.get(i).Y;
					System.out.printf( "%.5f %.5f %.5f\n", X1, X2, Y);
					
					double[] uAvAuBvBuCvC = new double[6];
					uAvAuBvBuCvC = FinduAvAuBvBuCvC();
					
					double[] output = new double[3];
					output = DerivativesOfOutputlayer(uAvAuBvBuCvC[5]);
										
					double[] HiddenLayer = new double[3];
					HiddenLayer = DerivativesOfHiddenLayer(uAvAuBvBuCvC[5], uAvAuBvBuCvC[0], 
							uAvAuBvBuCvC[2]);
									
					double[] derivativesOfWeights = DerivativesOfWeights(HiddenLayer, output, 
							uAvAuBvBuCvC[1], uAvAuBvBuCvC[3]);
					
					updateWEIGHTSsgd(derivativesOfWeights);
					for (int j = 0; j < WEIGHTS.length - 1; j++) {
						System.out.printf( "%.5f ", WEIGHTS[j]);
		            }
					System.out.printf( "%.5f", WEIGHTS[WEIGHTS.length - 1]);
					
					
					System.out.printf( "\n%.5f\n", setError(evaluation));				
				}
			    
			}
			else if (flag == 700) {
				int T = Integer.parseInt(args[11]);
				eta = Double.parseDouble(args[10]);
				WEIGHTS = new double[9];
				for (int i = 0; i < WEIGHTS.length; i++) {
	                WEIGHTS[i] = Double.parseDouble(args[i + 1]);
	            }
				ArrayList<setItem> training = FileReadSetItems(TRAINING);
				ArrayList<setItem> evaluation = FileReadSetItems(EVALUATION);
				
				for (int i = 0; i < T; i++) {
                    for (int j = 0; j < training.size(); j++) {
                    	X1 = training.get(j).X1;
    					X2 = training.get(j).X2;
    					Y = training.get(j).Y;
    					double[] uAvAuBvBuCvC = new double[6];
    					uAvAuBvBuCvC = FinduAvAuBvBuCvC();
    					
    					double[] output = new double[3];
    					output = DerivativesOfOutputlayer(uAvAuBvBuCvC[5]);
    										
    					double[] HiddenLayer = new double[3];
    					HiddenLayer = DerivativesOfHiddenLayer(uAvAuBvBuCvC[5], uAvAuBvBuCvC[0], 
    							uAvAuBvBuCvC[2]);
    									
    					double[] derivativesOfWeights = DerivativesOfWeights(HiddenLayer, output, 
    							uAvAuBvBuCvC[1], uAvAuBvBuCvC[3]);
    					
    					updateWEIGHTSsgd(derivativesOfWeights);
                    }
                    for (int j = 0; j < WEIGHTS.length - 1; j++) {
						System.out.printf( "%.5f ", WEIGHTS[j]);
		            }
					System.out.printf( "%.5f", WEIGHTS[WEIGHTS.length - 1]);
					System.out.printf( "\n%.5f\n", setError(evaluation));                  
                }
			}
			else if (flag == 800) {
				int T = Integer.parseInt(args[11]);
				eta = Double.parseDouble(args[10]);
				WEIGHTS = new double[9];
				for (int i = 0; i < WEIGHTS.length; i++) {
	                WEIGHTS[i] = Double.parseDouble(args[i + 1]);
	            }
				ArrayList<setItem> training = FileReadSetItems(TRAINING);
				ArrayList<setItem> evaluation = FileReadSetItems(EVALUATION);
				ArrayList<setItem> testdata = FileReadSetItems(TESTING);
				double previousSetError = Double.POSITIVE_INFINITY;
				int iteration = 0;
				
				for (int i = 0; i < T; i++) {
                    for (int j = 0; j < training.size(); j++) {
                    	X1 = training.get(j).X1;
    					X2 = training.get(j).X2;
    					Y = training.get(j).Y;
    					double[] uAvAuBvBuCvC = new double[6];
    					uAvAuBvBuCvC = FinduAvAuBvBuCvC();
    					
    					double[] output = new double[3];
    					output = DerivativesOfOutputlayer(uAvAuBvBuCvC[5]);
    										
    					double[] HiddenLayer = new double[3];
    					HiddenLayer = DerivativesOfHiddenLayer(uAvAuBvBuCvC[5], uAvAuBvBuCvC[0], 
    							uAvAuBvBuCvC[2]);
    									
    					double[] derivativesOfWeights = DerivativesOfWeights(HiddenLayer, output, 
    							uAvAuBvBuCvC[1], uAvAuBvBuCvC[3]);
    					
    					updateWEIGHTSsgd(derivativesOfWeights);
                    }
                   double setError = setError(evaluation);   
                    if (setError > previousSetError) {
                    	break;
                    } else {
                    	previousSetError = setError;
                    }
                    iteration++;
                }
				//how do i keep track of iteration?
				if (iteration == T) {
				System.out.println(iteration);
				} else {
					System.out.println(iteration + 1);
				}
				
				for (int j = 0; j < WEIGHTS.length - 1; j++) {
					System.out.printf( "%.5f ", WEIGHTS[j]);
	            }
				System.out.printf( "%.5f", WEIGHTS[WEIGHTS.length - 1]);
				System.out.printf( "\n%.5f\n", setError(evaluation));
				
				//calculating accuracy 
				double numOfCorrect = 0.0;
				for (int i = 0; i < testdata.size(); i++) {
					int testAnswer = (int) testdata.get(i).Y;
					int guess;
					
					X1 = testdata.get(i).X1;
					X2 = testdata.get(i).X2;
					
					if (FinduAvAuBvBuCvC()[5] >= 0.5) {
						guess = 1;
					} else {
						guess = 0;
					}
					
					if (testAnswer == guess) {
						numOfCorrect++;
					}
				}

				double accuracy = numOfCorrect / testdata.size();
				System.out.printf( "%.5f\n", accuracy);

			}
			
			
			else {
				throw new UnsupportedOperationException();
			}
			
			} catch (Exception e) {
		        System.out.print("$java Neural FLAG [arg1 arg2....arg12]");;
			}

	}

}
