package cubevb.transform;


public class LogExpM1 extends Transform {

	@Override
	public double forward(double x) {
		return softplus(x);
	}

	@Override
	public double backward(double z) {
        /* y = Log(Exp(x) - 1)
             = Log(1 - Exp(-x)) + x */
		double r = Math.log(1.0 - Math.exp(-z)) + z;
		if (parameter.getValue() != r) {
	        parameter.setValue(r);			
		}
        return r;
	}

	@Override
	public double jacobian_det(double z) {
		return -softplus(-z);
	}

	@Override
	public double backward(int dim, double oldValue, double newValue) {
		double r = Math.log(1.0 - Math.exp(-newValue)) + newValue;
		if (parameter.getValue(dim) != r) {
	        parameter.setValue(dim, r);			
		}
		return r;
	}

	@Override
	public double jacobian_det(int dim, double z) {
		return -softplus(-z);
	}

}
