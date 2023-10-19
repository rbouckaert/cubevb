package cubevb.transform;

public class LowerBound extends Transform {
    /* Transform from real line interval [a,inf] to whole real line. */
	double a;
    
    public LowerBound() {
    }
    
    public LowerBound(double a) {
    	this.a = a;
    }
    
    @Override
	public void initAndValidate() {
    	super.initAndValidate();
		this.a = parameter.getLower();
	}

    @Override
    public double backward(double z) {
    	double r = Math.exp(z) + a;
		if (parameter.getValue() != r) {
	        parameter.setValue(r);			
		}
		return r;
    }
    
    @Override
    public double forward(double x) {
    	return Math.log(x - a);
    }

    @Override
    public double jacobian_det(double x) {
    	return x;
    }

	@Override
	public double backward(int dim, double oldValue, double newValue) {
    	double r = Math.exp(newValue) + a;
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
