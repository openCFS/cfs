/**
 * this trigger update function is invoked by the "Update data" button and
 * start the 3 ajax request functions to request the result data, the viewed
 * results on the right and the iterative data
 */
function trigger_update() {
	if (!running) {
	
		running = true;
		update_results();
		update_view();
		update_iterations();
	} else {
		console.log("already running!");
	}
}

var running = false;

$( document ).ready(function() {
    $("#iteration * input").change(function () {
        console.log('input boxes changed');
        $("#iteration_num").html('-1');
    });
});

function update_results() {
	/*
	key = $("#simulation_id").html()

	var x = $("input[name=iteration_selector_x]:checked").val()

	var y = $("input[name=iteration_selector_y]:checked").map(function() {
		return $(this).val();
	}).get();

	iteration_num = $("#iteration_num").html();

	console.log(x)
	console.log(y)
	console.log(iteration_num)

	$.ajax({
		url : "/plot/" + key,
		data : {
			"x" : x,
			"y" : y,
			"iteration_num" : iteration_num
		},
		success : function(result) {

			// todo: stop refresh when simulation reached end
			if (result.indexOf("no_new_data") !== -1) {
				console.log("no new data!")
			} else {
				$("#iteration_plot").html(result);
			}

			if ($("input[name=enable_autorefresh]").prop('checked') == true) {
				setTimeout(function() {
					update_view()
				}, 1000);
			}
		},
		traditional : true,
		timeout : 10000
		// sets timeout to 10 seconds
	});*/
}

function update_view() {
	/*
	key = $("#simulation_id").html()

	var x = $("input[name=iteration_selector_x]:checked").val()

	var y = $("input[name=iteration_selector_y]:checked").map(function() {
		return $(this).val();
	}).get();

	iteration_num = $("#iteration_num").html();

	console.log(x)
	console.log(y)
	console.log(iteration_num)

	$.ajax({
		url : "/plot/" + key,
		data : {
			"x" : x,
			"y" : y,
			"iteration_num" : iteration_num
		},
		success : function(result) {

			// todo: stop refresh when simulation reached end
			if (result.indexOf("no_new_data") !== -1) {
				console.log("no new data!")
			} else {
				$("#iteration_plot").html(result);
			}

			if ($("input[name=enable_autorefresh]").prop('checked') == true) {
				setTimeout(function() {
					update_view()
				}, 1000);
			}
		},
		traditional : true,
		timeout : 10000
		// sets timeout to 10 seconds
	});*/
}


function update_iterations() {
	key = $("#simulation_id").html()

	var x = $("input[name=iteration_selector_x]:checked").val()

	var y1_it = $("input[name=iteration_selector_y1]:checked").map(function() {
		return $(this).val();
	}).get();
	
	var y2_it = $("input[name=iteration_selector_y2]:checked").map(function() {
		return $(this).val();
	}).get();

	var y1_res = $("input[name=result_selector_y1]:checked").map(function() {
		return $(this).val();
	}).get();
	
	var y2_res = $("input[name=result_selector_y2]:checked").map(function() {
		return $(this).val();
	}).get();

	iteration_num = $("#iteration_num").html();

	$.ajax({
		url : "/plot/" + key,
		data : {
			"x" : x,
			"y1_it" : y1_it,
			"y2_it" : y2_it,
			"y1_res" : y1_res,
			"y2_res" : y2_res,
			"iteration_num" : iteration_num
		},
		success : function(result) {

			// todo: stop refresh when simulation reached end
			if (result.indexOf("no_new_data") !== -1) {
				console.log("no new data!")
			} else {
				$("#iteration_plot").html(result);
			}

			if ($("input[name=enable_autorefresh]").prop('checked') == true) {
				setTimeout(function() {
					update_iterations()
				}, 1000);
			} else {
				running = false;
			}
		},
		error : function( jqXHR, textStatus, errorThrown ) {
			running = false;
			
			alert('There was an error ' + textStatus + 'during data request:\n' + errorThrown)
		},
		traditional : true,
		timeout : 10000
		// sets timeout to 10 seconds
	});
}
