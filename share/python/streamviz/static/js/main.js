/**
 * this trigger update function is invoked by the "Update data" button and
 * start the 3 ajax request functions to request the result data, the viewed
 * results on the right and the iterative data
 */
function trigger_update(continuous = false) {
	update_values(continuous);
	update_iterations(continuous);
}

function catalyst_send() {
	ip = $("#catalyst_ip").val();
	port = $("#catalyst_port").val();
	
	key = $("#simulation_id").html()
	
	/*$.ajax({
		url : "/catalyst_send/" + key,
		data : {
			"key" : key,
			"port" : port,
			"ip" : ip
		},
		success : function(result) {
		},
		error : function( jqXHR, textStatus, errorThrown ) {
		},
		traditional : true,
		timeout : 1000 // timeout of 1 second
	});*/
}

$( document ).ready(function() {

    // check the default x iterator
    $("input[name=iteration_selector_x][value=number]").prop('checked', true);
    
    // check the default objective
    if ($("#objective").length >= 1)  {
    	objective = $("#objective").html();
    	$("input[name=iteration_selector_y1][value=" + objective + "]").prop('checked', true);
    }
    
	if ($("#status").html().length == 0) {
		//nothing to do here
		return;
	}
	
	$("#catalyst_send_button").click(catalyst_send());
	
	$("#disable_results").click(function() {
		$("input[name='result_selector_view']").prop('checked', false);
		$("input[name='result_selector_view_bloch']").prop('checked', false);
	    trigger_update(false);
		//$("input[name='result_selector_view']:checked").buttonset('refresh');
	});
	
	continuous = false; // default for finished or abortet functions
	
	if ($("#status").html() == "running") {
		continuous = true
	}/* else {*/
		//we only need change on inpout changes when autorefresh is enabled
	    $("#iteration * input").change(function () {
	        console.log('input boxes changed');
	        update_iterations(false);
	    });

	    $("#logscale_y1").change(function () {
	        $("#iteration_num").html('-1');
	        update_iterations(false);
	    });
	    $("#logscale_y2").change(function () {
	        $("#iteration_num").html('-1');
	        update_iterations(false);
	    });
	    $("input[name='result_selector_view']").change(function () {
	        $("#iteration_num").html('-1');
	        update_iterations(false);
	    });
	    $("input[name='result_selector_view_bloch']").change(function () {
	        $("#iteration_num").html('-1');
	        update_iterations(false);
	    });
	//}

    trigger_update(continuous);
});

function update_values(refresh = false) {
	key = $("#simulation_id").html()
	
	iteration_num = $("#iteration_num").html();

	if ($("#status").html() == "finished") {
		iteration_num = -1;
	}

	$.ajax({
		url : "/values/" + key,
		data : {
			"key" : key,
			"iteration_num" : iteration_num
		},
		success : function(result) {

			if (result.indexOf("no_new_data") !== -1) {
				console.log("no new data!")
			} else {
				// if new data, then also send catalyst
				if ($('#catalyst_send_auto_update').is(":checked")) {
					catalyst_send();
				}
				data = JSON.parse(result);

				for (var name in data['iterations']) {
					$('#td_seqit_container_' + name).html(data['iterations'][name]);
				}

				for (var name in data['results']) {
					$('#td_seqres_container_' + name).html(data['results'][name]);
				}
			}

			
			if (refresh) {
				if ($("#status").html() != "running") {
					window.location.reload(false); 
				}
				else {
					setTimeout(function() {
						update_values(true)
					}, 300);
				}
			}
		},
		error : function( jqXHR, textStatus, errorThrown ) {
			alert('There was an error ' + textStatus + 'during data request:\n' + errorThrown);
		},
		traditional : true,
		timeout : 10000
		// sets timeout to 10 seconds
	});
}

function update_iterations(refresh = false) {
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

	if ($("#status").html() == "finished") {
		iteration_num = -1;
	}
	
	if (!refresh) {
		iteration_num = -1;
	}

	logscale_y1 = false;
	if ($("#logscale_y1").prop('checked') == true) {
		logscale_y1 = true;
	}

	logscale_y2 = false;
	if ($("#logscale_y2").prop('checked') == true) {
		logscale_y2 = true;
	}

	show_bloch = 'false';
	$("input[name='result_selector_view_bloch']:checked").each( function () {
		show_bloch = 'true';
	});

	result_array = [];
	$("input[name='result_selector_view']:checked").each( function () {
		result_array.push($(this).val());
	});

	$.ajax({
		url : "/plot/" + key,
		data : {
			"x" : x,
			"y1_it" : y1_it,
			"y2_it" : y2_it,
			"y1_res" : y1_res,
			"y2_res" : y2_res,
			"iteration_num" : iteration_num,
			"logscale_y1" : logscale_y1,
			"logscale_y2" : logscale_y2,
			"view_results" : result_array,
			"view_result_bloch" : show_bloch
		},
		success : function(result) {

			// todo: stop refresh when simulation reached end
			if (result.indexOf("no_new_data") !== -1) {
				console.log("no new data!")
			} else {
				$("#iteration_plot").html(result);
			}

			if (refresh && !($("#status").html() != "running")) {
				setTimeout(function() {
					update_iterations(true)
				}, 300);
			}
		},
		error : function( jqXHR, textStatus, errorThrown ) {
			alert('There was an error ' + textStatus + 'during data request:\n' + errorThrown);
		},
		traditional : true,
		timeout : 10000
		// sets timeout to 10 seconds
	});
}
