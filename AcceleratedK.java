import org.jblas.DoubleMatrix;
import java.lang.Math;
import org.jblas.MatrixFunctions;
import org.jblas.ranges.IntervalRange;

import java.util.*;

public class AcceleratedK {
	public static DoubleMatrix calcDistances(DoubleMatrix data, DoubleMatrix center){
		int dataRows = data.getRows();
		int centerRows = center.getRows();
		int dataCols = data.getColumns();
		int centerCols = center.getColumns();
		DoubleMatrix distances = null;
		if(centerRows == 1){
			distances = MatrixFunctions.pow(data, 2).rowSums().sub(data.mul(center.transpose()).mmul(2)).add(
								center.mul(center.transpose()));		
		}
		else{
			//assert centerRows and dataRows same
			distances = MatrixFunctions.pow(data.sub(center), 2).rowSums();
		}
		return distances;
	}
	public static DoubleMatrix allDistances(DoubleMatrix centers){
		int k = centers.getRows();
		DoubleMatrix centDist = DoubleMatrix.zeros(k,k);
		for(int j = 1;j < k;j++){
			centDist.put(new IntervalRange(0,j), new IntervalRange(j,j+1),
						calcDistances(centers.get(new IntervalRange(0,j), new IntervalRange(0,k)),
										centers.getRow(j)));
		}
		centDist = centDist.add(centDist.transpose());
		return centDist;
	}
}
