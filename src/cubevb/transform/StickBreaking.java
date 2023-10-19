package cubevb.transform;


import java.util.Arrays;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;;

@Description("Transforms K - 1 dimensional simplex space (k values in [0,1] and that sum to 1) to a K - 1 vector of real values. ")
//    Primarily borrowed from the STAN implementation.
public class StickBreaking extends Transform {
	public final Input<Double> sumInput = new Input<>("sum", "total sum of individual components", 1.0);

	double [] shape;
	
	double sum;

	// A small value for numerical stability in invlogit.
    double eps;

    // dimension of parameter 
    int k, Km1;
    double [] eq_share;
    
    public StickBreaking() {
    }
    
	public StickBreaking(RealParameter realParameter) {
		initByName("parameter", realParameter);
	}

	public StickBreaking(RealParameter realParameter, double sum) {
		initByName("parameter", realParameter, "sum", sum);
	}

	@Override
	public void initAndValidate() {
		super.initAndValidate();
		
		sum = sumInput.get();
		
		k = parameter.getDimension();
        Km1 = k - 1;
        eq_share = new double[k-1];
        for (int i = 0; i < k-1; i++) {
    		// implement based on formulas at https://mc-stan.org/docs/2_22/reference-manual/simplex-transform-section.html
        	eq_share[i] = Math.log(1.0 / (Km1 - i));
        }
	}

	@Override
	public double[] forward(double[] x) {
        double [] s = new double[Km1];
        s[0] = 1.0-x[0]/sum;
        for (int i = 1; i < Km1; i++) {
        	s[i] = s[i-1] - x[i]/sum;
        }
        // z = x0 / s
        double [] z = new double[Km1];
        z[0] = 1.0-s[0];
        for (int i = 1; i < s.length; i++) {
        	z[i] = 1-s[i] / s[i-1];
        }
        
        double [] y = new double[Km1];
        for (int i = 0; i < Km1; i++) {
        	y[i] = Transform.logit(z[i]) - eq_share[i];
        }
        return y;
	}
	
	@Override
	public double[] backward(double[] y) {
        double [] z = new double[k-1];
        for (int i = 0; i < Km1; i++) {
            z[i]= Transform.invlogit(y[i] + eq_share[i]);
        }
        
        double [] x = new double[k];
        
        double s = 0;
        for (int i = 0; i < Km1; i++) {
        	x[i] = (sum - s) * z[i];
        	s += x[i];
    		if (parameter.getValue(i) != x[i]) {
    			parameter.setValue(i, x[i]);
    		}
        }
        x[Km1] = sum - s;
		if (parameter.getValue(Km1) != x[Km1]) {
			parameter.setValue(Km1, x[Km1]);
		}
        
        return x;	
	}

	@Override
	public double jacobian_det(double[] y) {
		// implement based on formulas at https://mc-stan.org/docs/2_22/reference-manual/simplex-transform-section.html
        double [] z = new double[k-1];
        for (int i = 0; i < Km1; i++) {
            z[i]= Transform.invlogit(y[i] + eq_share[i]);
        }
        
        double [] x = new double[k];
        double s = 0;
        for (int i = 0; i < Km1; i++) {
        	x[i] = (sum - s) * z[i];
        	s += x[i];
        }
        x[Km1] = 1 - s;

        double [] S = new double[k];
        S[0] = 0;
        for (int i = 0; i < Km1; i++) {
        	S[i+1] = S[i] + x[i]/sum;
        }
        
        double sum = 0;
        for (int i = 0; i < Km1; i++) {
        	sum += Math.log1p(-S[i]);
            sum += Math.log(z[i]);
            sum += Math.log1p(-z[i]);
        }
        return sum;
	}
	
	
	@Override
	public int getDimension() {
		return parameter.getDimension() - 1;
	} 

	
	public static void main(String[] args) {
		RealParameter p = new RealParameter("0.2 0.2 0.2 0.4");
		StickBreaking b = new StickBreaking(p);
		double [] y = new double[]{0.0, 0.0, 0.0};
		double [] x = b.backward(y);
		double [] y2 = b.forward(x);
		
		
		
		System.out.println("y="+Arrays.toString(y));
		System.out.println("y2="+Arrays.toString(y2));
		

		
		// test how backward(forward(x)) differs from x
		// should be O(1e-16)
		double maxDiff = 0;
		Randomizer.setSeed(1237);
		for (int i = 0; i < 10000; i++) {
			double [] values = new double[4];
			for (int j = 0; j < values.length - 1; j++) {
				values[j] = Randomizer.nextDouble();
			}
			values[values.length-1] = 1.0;
			Arrays.sort(values);
			for (int j = values.length - 1; j > 0; j--) {
				values[j] = values[j] - values[j-1];
			}			
			
			y2 = b.forward(values);
			x = b.backward(y2);
			for (int j = 0; j < values.length; j++) {
				if (Math.abs(x[j] - values[j]) > maxDiff) {
					maxDiff = Math.abs(x[j] - values[j]);
					System.err.println(i + " " + maxDiff);
				}
			}
		}
		
		
		System.out.println("x="+Arrays.toString(x));
		System.out.println(Arrays.toString(p.getDoubleValues()));
		System.out.println("J="+b.jacobian_det(y));

		p = new RealParameter("0.4 0.4 0.4 0.8");
		b = new StickBreaking(p, 2.0);
		y = new double[]{0.1, 0.1, 0.0};
		x = b.backward(y);
		y2 = b.forward(x);
		
		System.out.println("y="+Arrays.toString(y));
		System.out.println("y2="+Arrays.toString(y2));
		
		System.out.println("x="+Arrays.toString(x));
		System.out.println(Arrays.toString(p.getDoubleValues()));
		System.out.println("J="+b.jacobian_det(y));
		
		double sum = 0;
		for (double d : p.getValues()) {
			sum += d;
		}
		System.out.println("sum = " + sum);
	
	}

	@Override
	public double forward(double x) {
		throw new RuntimeException("Cannot use single dimension call on StickBreaking transform");
	}

	@Override
	public double backward(double z) {
		throw new RuntimeException("Cannot use single dimension call on StickBreaking transform");
	}

	@Override
	public double jacobian_det(double x) {
		throw new RuntimeException("Cannot use single dimension call on StickBreaking transform");
	}

	@Override
	public double backward(int dim, double oldValue, double newValue) {
		double [] v = forward();
		v[dim] = newValue;
		double [] x = backward(v);
		return x[dim];
	}

	@Override
	public double jacobian_det(int dim, double x) {
		double [] v = forward();
		v[dim] = Double.NEGATIVE_INFINITY;
		double j = jacobian_det(v);
		v[dim] = x;
		j = jacobian_det(v) - j;
		return j;
	}
	
}

