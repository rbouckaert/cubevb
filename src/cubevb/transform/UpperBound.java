package cubevb.transform;

public class UpperBound extends Transform {
    /* Transform from real line interval [inf,b] to whole real line. */
	double b;
    
    public UpperBound() {    	
    }
    
    public UpperBound(double b) {
    	this.b = b;
    }
    
    
    @Override
	public void initAndValidate() {
    	super.initAndValidate();
		this.b = parameter.getUpper();
	}


    @Override
    public double backward(double z) {
    	double r = b - Math.exp(z);
		if (parameter.getValue() != r) {
	        parameter.setValue(r);			
		}
		return r;
    }
    
    @Override
    public double forward(double x) {
    	return Math.log(b - x);
    }

    @Override
    public double jacobian_det(double x) {
    	return x;
    }

	@Override
	public double backward(int dim, double oldValue, double newValue) {
    	double r = b - Math.exp(newValue);
		if (parameter.getValue(dim) != r) {
	        parameter.setValue(dim, r);			
		}
		return r;
	}

	@Override
	public double jacobian_det(int dim, double x) {
		return x;
	}
}	
