function openClose() {
	if ($("#post-box").css("display") == "block") {
		$("#post-box").hide();
		$("#btn-post-box").text("Open input box");
	} else {
		$("#post-box").show();
		$("#btn-post-box").text("Close input box");
	}
}



function postItem() {
	let gene_seq = $('#gene_seq').val();
	let oligomer_size = $('#oligomer_size').val();
	let overlap_size = $('#overlap_size').val();
	let optimal_temp = $('#optimal_temp').val();
	let temp_range = $('#temp_range').val();
	let cluster_size = $('#cluster_size').val();
	let cluster_range = $('#cluster_range').val();

	$.ajax({
		type: "POST",
		url: "/insert_item",
		data: { 
			gene_seq : gene_seq,
			oligomer_size : oligomer_size,
			overlap_size : overlap_size,
			optimal_temp : optimal_temp,
			temp_range : temp_range,
			cluster_size : cluster_size,
			cluster_range : cluster_range
		},
		success: function (response) { 
			if (response["result"] == "success") {
				//alert(response["msg"]);
				window.location.reload();
				//PrettyFasta()
			}
		}
	})
}

function old_listItems() {
	$.ajax({
		type: "GET",
		url: "/list_items",
		data: {},
		success: function (response) {
			if (response["result"] == "success") {
				let items = response['items'];
				//console.log(items);
				for (let i = 0; i < items.length; i++){
					let gene_seq = items[i]['gene_seq'];
					let oligomer_size = items[i]['oligomer_size'];
					let overlap_size = items[i]['overlap_size']; 
					let optimal_temp = items[i]['optimal_temp'];
					let temp_html = `
						<div>
						<p class="fasta">${gene_seq}</p>
						</div>
						`;
					$('#item_list').append(temp_html);
				}
			}
		}
	})
}


function old2_listItems(){
	$(document).ready(function () {
		 $.noConflict();
		$("#item_list").html("");

		console.log("test")
		$('#item-list-table').DataTable({
			"ajax": {
				"url": "/list_items",
				"datatype": "json",
				"type": "GET",
				"dataSrc": "data",
				"contentType": "application/json"
			},
			"columns": [
				{ "data": "gene_seq_short" }, 
				{ "data": "oligomer_size" }, 
				{ "data": "optimal_temp" }, 
				{ "data": "cluster_size" } 
			],
			"select": { 
				style: 'single' 
			}
		});
	});
}

function detailFormat(d) {
	var s = "";
	s += "<div class='target-sequence-label'>Target full sequence</div>";
	s += "<div class='fasta'>"+d.gene_seq+"</div>";
	s += "<button id='mybutton' type='button' class='btn btn-link' data-toggle='modal' data-target='#exampleModalCenter'> Show design results </button>"
	return s
}

function listItems(){
	$(document).ready(function() {
		$.noConflict();
		var dt = $('#item-list-table').DataTable({
			"ajax": {
				"url": "/list_items",
				"datatype": "json",
				"type": "GET",
				"dataSrc": "data",
				"contentType": "application/json"
			},
			"columns": [
				{
					"class": "details-control",
					"orderable": false,
					"data": null,
					"defaultContent": ""
				},
				{ "data": "gene_seq_short" }, 
				{ "data": "oligomer_size" }, 
				{ "data": "overlap_size" }, 
				{ "data": "optimal_temp" }, 
				{ "data": "cluster_size" } 
			],
			"select": { 
				style: 'single' 
			},
			"order": [[0, 'asc']]
		});

		var detailRows = [];
		$('#item-list-table tbody').on( 'click', 'tr td.details-control', function () {
			var tr = $(this).closest('tr');
			var row = dt.row( tr );
			var idx = $.inArray( tr.attr('id'), detailRows );

			if ( row.child.isShown() ) {
				tr.removeClass( 'details' );
				row.child.hide();
				detailRows.splice( idx, 1 );
			}else{
				tr.addClass('details');
				row.child(detailFormat(row.data())).show();
				PrettyFasta();
				if(idx == 1){
					detailRows.push( tr.attr('id') );
				}
			}
		});

		dt.on( 'draw', function(){
			$.each( detailRows, function(i, id) {
				$('#'+id+' td.details-control').trigger( 'click' );
			});
		});
	});
}



function showItems(){
	var myModal = document.getElementById('myModal')
	var myInput = document.getElementById('myInput')

	myModal.addEventListener('shown.bs.modal', function () {
		  myInput.focus()
	})
}

$( window ).on("load", listItems())








