package cubevb.transform;

public class Tuple implements Comparable<Tuple> {
	int i, j;
	double covar;
	
	public Tuple(int i, int j, double covar) {
		this.i = i;
		this.j = j;
		this.covar = covar;
	}
	
	@Override
	public int compareTo(Tuple o) {
		if (o.covar > covar) {
			return 1;
		} else if (o.covar < covar) {
			return -1;
		}
		return 0;
	}
}
