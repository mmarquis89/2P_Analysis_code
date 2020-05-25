
test = EventAlignedData(expList(1:40, :));
test = test.set_align_event('odor');
test = test.extract_roi_data();
test = test.extract_fictrac_data();

% test = test.create_filter_event_vectors();
test2 = test.output_analysis_subset([4 4], ');

% head(test.sourceMd, 1)
% head(test.alignEventTable, 1)
% head(test.eventFlTable, 1)
% head(test.eventFilterTable, 1)
