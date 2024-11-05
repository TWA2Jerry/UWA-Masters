using CairoMakie
scaling_factor = 2.5

arrow_marker = BezierPath([
       MoveTo(Point(0,0)),
       LineTo(Point(-0.3, -0.3)*scaling_factor),
           LineTo(Point(-0.3, -0.10)*scaling_factor),
	   LineTo(Point(-1.5, -0.1)*scaling_factor),
	   LineTo(Point(-1.5, 0.1)*scaling_factor),
		LineTo(Point(-0.3, 0.1)*scaling_factor),
       LineTo(Point(-0.3, 0.3)*scaling_factor),
       ClosePath()
       ])

arrow_marker_curved = BezierPath([
       MoveTo(Point(0,0)),
       LineTo(Point(-0.3, -0.3)*scaling_factor),
           LineTo(Point(-0.3, -0.10)*scaling_factor),
       CurveTo(Point(-0.65, -0.225)*scaling_factor, Point(-0.65, -0.225)*scaling_factor, Point(-1.5, 0.0)*scaling_factor),
       CurveTo(Point(-0.65, 0.225)*scaling_factor, Point(-0.65, 0.225)*scaling_factor, Point(-0.3, 0.10)*scaling_factor),
       LineTo(Point(-0.3, 0.3)*scaling_factor),
       ClosePath()
       ])
