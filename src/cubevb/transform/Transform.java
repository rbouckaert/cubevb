package cubevb.transform;

import java.util.*;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.CompoundRealParameter;
import beast.base.inference.parameter.RealParameter;

// based on pymc implementation from 
// https://github.com/pymc-devs/pymc3/blob/master/pymc3/distributions/transforms.py

abstract public class Transform extends BEASTObject {
//	public final Input<List<RealParameter>> parameterInput = new Input<>("parameter", "parameter that is the output of stick breaking construction", new ArrayList<>(), Validate.REQUIRED);
	public final Input<RealParameter> parameterInput = new Input<>("parameter", "parameter that is the output of stick breaking construction", Validate.REQUIRED);
	
    RealParameter parameter;
    
    public Transform() {}
    public Transform(RealParameter parameter) {initByName("parameter", parameter);}
	
    @Override
	public void initAndValidate() {
    	parameter = parameterInput.get();
//    	if (parameterInput.get().size() == 1) {
//    		parameter = parameterInput.get().get(0);
//    	} else {
//    		parameter = new CompoundRealParameter();
//    		parameter.initByName("parameter", parameterInput.get());
//    		String id = "";
//    		for (RealParameter p : parameterInput.get()) {
//    			id += p.getID() + "-";
//    		}
//    		id = id.substring(0, id.length() - 1);
//    		parameter.setID(id);
//    	}
	}

    public RealParameter getParameter() {
    	return parameter;    	
    }
    
	/* A transformation of a random variable from one space into another.
    Attributes
    ----------
    name : str
    */ 

    // name = ""

	public int getDimension() {
		return parameter.getDimension();
	} 
	
    public double [] forward() {return forward(parameter.getDoubleValues());}
    public double [] forward(double [] x) {
    	double [] y = new double[x.length];
    	for (int i = 0; i < x.length; i++) {
    		y[i] = forward(x[i]);
    	}
    	return y;
    }
    abstract public double forward(double x);
        /* Applies transformation forward to input variable `x`.
        When transform is used on some distribution `p`, it will transform the random variable `x` after sampling
        from `p`.
        Parameters
        ----------
        x : tensor
            Input tensor to be transformed.
        Returns
        --------
        tensor
            Transformed tensor.
        */ 


    public double [] backward(double [] z) {
    	double [] x = new double[z.length];
    	for (int i = 0; i < z.length; i++) {
    		x[i] = backward(i, Double.NaN, z[i]);
    	}
    	for (int i = 0; i < z.length; i++) {
    		parameter.setValue(i, x[i]);
    	}
    	return x;
    }
    abstract public double backward(double z);
        /* Applies inverse of transformation to input variable `z`.
        When transform is used on some distribution `p`, which has observed values `z`, it is used to
        transform the values of `z` correctly to the support of `p`.
        Parameters
        ----------
        z : tensor
            Input tensor to be inverse transformed.
        Returns
        -------
        tensor
            Inverse transformed tensor.
        */ 

    abstract public double backward(int dim, double oldValue, double newValue);

    public double jacobian_det(double [] x) {
    	double jac = 0;
    	for (double d : x) {
    		jac += jacobian_det(d);
    	}
    	return jac;
    }
    abstract  public double jacobian_det(double x);
        /* Calculates logarithm of the absolute value of the Jacobian determinant for input `x`.
        Parameters
        ----------
        x : tensor
            Input to calculate Jacobian determinant of.
        Returns
        -------
        tensor
            The log abs Jacobian determinant of `x` w.r.t. this transform.
        */
    
    abstract  public double jacobian_det(int dim, double x);
 	
	static public double softplus(double x) {
		final double z = Math.exp(x);
		final double y = 1.0 + z;
		if (y > 1.0) {
			return Math.log(y);
		}
		return z;
	}
	
	static public double softplusInverse(double x) {
		final double y = Math.exp(x) - 1.0;
		if (y > 0) {
			return Math.log(y);
		}
		return Math.log(x);
	}
	
	static public double logit(double x) {
		return Math.log(x/(1.0-x));
	}
	
	static public double invlogit(double x) {
		return 1.0/(1+Math.exp(-x));
	}

	static public double sigmoid(double x) { // = logit(x)
		return 1.0/(1+Math.exp(-x));
	}
	
	/** node that is transformed **/
	public StateNode getTransformNode() {
		return parameter;
	}
	
	/** node that is in the state **/
	public StateNode getStateNode() {
		return parameter;
	}
}

