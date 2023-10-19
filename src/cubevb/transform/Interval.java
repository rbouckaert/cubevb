package cubevb.transform;

public class Interval extends Transform {
    /*Transform from real line interval [a,b] to whole real line.*/
	
	double a, b;
	
	public Interval() {
	}
	
	public Interval(double a, double b) {
		this.a = a;
		this.b = b;
	}
	
	
    @Override
	public void initAndValidate() {
    	super.initAndValidate();
		this.a = parameter.getLower();
		this.b = parameter.getUpper();
	}
	
    @Override
    public double backward(double z) {
        double r = (b - a) * sigmoid(z) + a;
		if (parameter.getValue() != r) {
	        parameter.setValue(r);			
		}
        return r;
    }

    @Override
    public double forward(double x) {
    	return Math.log(x - a) - Math.log(b - x);
    }
    
    @Override
    public double jacobian_det(double x) {
    	double s = softplus(-x);
    	return Math.log(b - a) - 2 * s - x;
    }

	@Override
	public double backward(int dim, double oldValue, double x) {
        double r = (b - a) * sigmoid(x) + a;
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
