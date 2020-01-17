//use Random;

proc main(args: [] string){
	var count: int = 0;//sync int = 0; // number of points in the 1st quadrant of unit circle
	var N: int = 1000;
	if args.size >= 1 {
		N = args(1): int;
	}

	//var randStream = new owned RandomStream(real, 0);
	for i in 0..N {
		var x: real = 0.0;//randStream.getNext();
		var y: real = 0.0;//randStream.getNext();
		var z: real = x * x + y * y;
		writeln(x,y,z);
		if z <= 1 {
			count += 1;
		}
	}

	var pi: real = count / N * 4;
	writeln("calculated value of pi:",pi);
}
